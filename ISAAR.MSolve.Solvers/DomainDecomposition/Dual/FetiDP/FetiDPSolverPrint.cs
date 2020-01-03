using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Rigid body modes do not have to be computed each time the stiffness matrix changes. 
//TODO: Optimizations for the case that stiffness changes, but connectivity remains the same!
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP
{
    public class FetiDPSolverPrint : ISolver
    {
        internal const string name = "FETI-DP Solver"; // for error messages
        private readonly IFetiDPCoarseProblemSolver coarseProblemSolver;
        private readonly ICornerNodeSelection cornerNodeSelection;
        private readonly ICrosspointStrategy crosspointStrategy = new FullyRedundantConstraints();
        private readonly IDofOrderer dofOrderer;
        private readonly FetiDPDofSeparator dofSeparator;
        private readonly IFetiDPInterfaceProblemSolver interfaceProblemSolver;
        private readonly Dictionary<int, IFetiDPSubdomainMatrixManager> matrixManagers;
        private readonly Dictionary<int, IFetiSubdomainMatrixManager> matrixManagersGeneral; //TODO: redesign. They are the same as above, but Dictionary is not covariant
        private readonly Dictionary<int, ISingleSubdomainLinearSystem> linearSystems;
        private readonly IStructuralModel model;
        private readonly IFetiPreconditionerFactory preconditionerFactory;
        private readonly bool problemIsHomogeneous;
        private readonly IStiffnessDistribution stiffnessDistribution;

        //TODO: fix the mess of Dictionary<int, ISubdomain>, List<ISubdomain>, Dictionary<int, Subdomain>, List<Subdomain>
        //      The concrete are useful for the preprocessor mostly, while analyzers, solvers need the interface versions.
        //      Lists are better in analyzers and solvers. Dictionaries may be better in the preprocessors.
        private readonly Dictionary<int, ISubdomain> subdomains;

        private bool factorizeInPlace = true;
        private FetiDPFlexibilityMatrix flexibility;
        private bool isStiffnessModified = true;
        private FetiDPLagrangeMultipliersEnumerator lagrangeEnumerator;
        private IFetiPreconditioner preconditioner;
        private FetiDPSubdomainGlobalMapping subdomainGlobalMapping;

        private FetiDPSolverPrint(IStructuralModel model, ICornerNodeSelection cornerNodeSelection,
            IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory, IDofOrderer dofOrderer, 
            IFetiPreconditionerFactory preconditionerFactory, bool problemIsHomogeneous, 
            IFetiDPInterfaceProblemSolver interfaceProblemSolver)
        {
            // Model
            if (model.Subdomains.Count == 1) throw new InvalidSolverException(
                $"{name} cannot be used if there is only 1 subdomain");
            this.model = model;
            this.cornerNodeSelection = cornerNodeSelection;

            // Subdomains
            subdomains = new Dictionary<int, ISubdomain>();
            foreach (ISubdomain subdomain in model.Subdomains) subdomains[subdomain.ID] = subdomain;

            // Matrix managers and linear systems
            matrixManagers = new Dictionary<int, IFetiDPSubdomainMatrixManager>();
            matrixManagersGeneral = new Dictionary<int, IFetiSubdomainMatrixManager>();
            this.linearSystems = new Dictionary<int, ISingleSubdomainLinearSystem>();
            var externalLinearSystems = new Dictionary<int, ILinearSystem>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                var matrixManager = matrixManagerFactory.CreateMatricesManager(subdomain);
                matrixManagers[s] = matrixManager;
                matrixManagersGeneral[s] = matrixManager;
                this.linearSystems[s] = matrixManager.LinearSystem;
                externalLinearSystems[s] = matrixManager.LinearSystem;

                //TODO: This will call HandleMatrixWillBeSet() once for each subdomain. For now I will clear the data when 
                //      BuildMatrices() is called. Redesign this.
                //matrixManager.LinearSystem.MatrixObservers.Add(this); 
            }
            LinearSystems = externalLinearSystems;

            this.dofOrderer = dofOrderer;
            this.dofSeparator = new FetiDPDofSeparator();
            this.preconditionerFactory = preconditionerFactory;

            // Interface problem
            this.coarseProblemSolver = matrixManagerFactory.CreateCoarseProblemSolver(model.Subdomains);
            this.interfaceProblemSolver = interfaceProblemSolver;

            // Homogeneous/heterogeneous problems
            this.problemIsHomogeneous = problemIsHomogeneous;
            if (problemIsHomogeneous) this.stiffnessDistribution = new HomogeneousStiffnessDistribution(model, dofSeparator);
            else this.stiffnessDistribution = new HeterogeneousStiffnessDistribution(model, dofSeparator);
        }

        public Dictionary<int, HashSet<INode>> CornerNodesOfSubdomains { get; private set; }
        public IReadOnlyDictionary<int, ILinearSystem> LinearSystems { get; }
        public SolverLogger Logger { get; } = new SolverLogger(name);
        public string Name => name;

        public Dictionary<int, IMatrix> BuildGlobalMatrices(IElementMatrixProvider elementMatrixProvider)
        {
            HandleMatrixWillBeSet(); //TODO: temporary solution to avoid this getting called once for each linear system/observable

            var watch = new Stopwatch();
            watch.Start();
            var matrices = new Dictionary<int, IMatrix>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                int s = subdomain.ID;
                IMatrix Kff;
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine($"{this.GetType().Name}: Assembling the free-free stiffness matrix of subdomain {s}");
                    Kff = matrixManagers[s].BuildGlobalMatrix(subdomain.FreeDofOrdering,
                        subdomain.Elements, elementMatrixProvider);
                    linearSystems[s].Matrix = Kff; //TODO: This should be done by the solver not the analyzer. This method should return void.
                }
                else
                {
                    Kff = (IMatrix)(linearSystems[s].Matrix); //TODO: remove the cast
                }
                matrices[s] = Kff;
            }
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);

            this.Initialize(); //TODO: Should this be called by the analyzer? Probably not, since it must be called before DistributeBoundaryLoads().
            return matrices;
        }

        //TODO: There is a lot of repetition between this method and BuildGlobalMatrices
        public Dictionary<int, (IMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr)> BuildGlobalSubmatrices(IElementMatrixProvider elementMatrixProvider)
        {
            HandleMatrixWillBeSet(); //TODO: temporary solution to avoid this getting called once for each linear system/observable

            var watch = new Stopwatch();
            watch.Start();
            var matrices = new Dictionary<int, (IMatrix Aff, IMatrixView Afc, IMatrixView Acf, IMatrixView Acc)>();
            foreach (ISubdomain subdomain in model.Subdomains) //TODO: this must be done in parallel
            {
                int s = subdomain.ID;
                if (!subdomain.StiffnessModified)
                {
                    throw new NotImplementedException("This optimization is not implemented");
                }
                if (subdomain.ConstrainedDofOrdering == null)
                {
                    throw new InvalidOperationException("In order to build the matrices corresponding to constrained dofs of,"
                        + $" subdomain {s}, they must have been ordered first.");
                }
                (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) =
                    matrixManagers[s].BuildGlobalSubmatrices(subdomain.FreeDofOrdering, subdomain.ConstrainedDofOrdering, 
                    subdomain.Elements, elementMatrixProvider);
                matrices[s] = (Kff, Kfc, Kcf, Kcc);
                linearSystems[s].Matrix = Kff; //TODO: This should be done by the solver not the analyzer. This method should return void.
            }
            watch.Stop();
            Logger.LogTaskDuration("Matrix assembly", watch.ElapsedMilliseconds);

            this.Initialize(); //TODO: Should this be called by the analyzer? Probably not, since it must be called before DistributeBoundaryLoads().
            return matrices;
        }

        //TODO: this and the fields should be handled by a class that handles dof mappings.
        public Dictionary<int, SparseVector> DistributeNodalLoads(Table<INode, IDofType, double> globalNodalLoads)
            => subdomainGlobalMapping.DistributeNodalLoads(subdomains, globalNodalLoads);

        //TODO: this and the fields should be handled by a class that handles dof mappings.
        public Vector GatherGlobalDisplacements(Dictionary<int, IVectorView> subdomainDisplacements)
            => subdomainGlobalMapping.GatherGlobalDisplacements(subdomainDisplacements);

        public void HandleMatrixWillBeSet()
        {
            isStiffnessModified = true;
            foreach (ISubdomain subdomain in subdomains.Values)
            {
                if (subdomain.StiffnessModified)
                {
                    Debug.WriteLine($"{this.GetType().Name}: Clearing saved matrices of subdomain {subdomain.ID}.");
                    matrixManagers[subdomain.ID].Clear();
                }
            }
            flexibility = null;
            preconditioner = null;
            coarseProblemSolver.ClearCoarseProblemMatrix();

            //stiffnessDistribution = null; //WARNING: do not dispose of this. It is updated when BuildGlobalMatrix() is called.
        }

        public void Initialize()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Identify corner nodes
            CornerNodesOfSubdomains = cornerNodeSelection.SelectCornerNodesOfSubdomains(); //TODO: Could this cause change in connectivity?

            // Define boundary / internal dofs
            dofSeparator.DefineGlobalBoundaryDofs(model, CornerNodesOfSubdomains);
            dofSeparator.DefineGlobalCornerDofs(model, CornerNodesOfSubdomains);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                HashSet<INode> cornerNodes = CornerNodesOfSubdomains[s];
                if (subdomain.ConnectivityModified)
                {
                    //TODO: should I cache this somewhere?
                    //TODO: should I use subdomain.Nodes.Except(cornerNodes) instead?
                    IEnumerable<INode> remainderAndConstrainedNodes = subdomain.Nodes.Where(node => !cornerNodes.Contains(node));

                    Debug.WriteLine($"{this.GetType().Name}: Separating and ordering corner-remainder dofs of subdomain {s}");
                    dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes, remainderAndConstrainedNodes);

                    Debug.WriteLine($"{this.GetType().Name}: Reordering internal dofs of subdomain {s}.");
                    matrixManagers[s].ReorderRemainderDofs(dofSeparator, subdomain);

                    Debug.WriteLine($"{this.GetType().Name}: Separating and ordering boundary-internal dofs of subdomain {s}");
                    dofSeparator.SeparateBoundaryInternalDofs(subdomain, remainderAndConstrainedNodes);
                }
            }

            // This must be called after determining the corner dofs of each subdomain
            //TODO: Perhaps the corner dofs of each subdomain should be handled in dofSeparator.DefineGlobalCornerDofs();
            coarseProblemSolver.ReorderCornerDofs(dofSeparator);
            dofSeparator.CalcCornerMappingMatrices(model, CornerNodesOfSubdomains);

            //TODO: B matrices could also be reused in some cases
            // Define lagrange multipliers and boolean matrices. 
            this.lagrangeEnumerator = new FetiDPLagrangeMultipliersEnumerator(crosspointStrategy, dofSeparator);
            if (problemIsHomogeneous) lagrangeEnumerator.DefineBooleanMatrices(model); // optimization in this case
            else lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model);

            // Log dof statistics
            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            int numExpandedDomainFreeDofs = 0;
            foreach (var subdomain in model.Subdomains)
            {
                numExpandedDomainFreeDofs += subdomain.FreeDofOrdering.NumFreeDofs;
            }
            Logger.LogNumDofs("Expanded domain dofs", numExpandedDomainFreeDofs);
            Logger.LogNumDofs("Lagrange multipliers", lagrangeEnumerator.NumLagrangeMultipliers);
            Logger.LogNumDofs("Corner dofs", dofSeparator.NumGlobalCornerDofs);

            // Use the newly created stiffnesses to determine the stiffness distribution between subdomains.
            //TODO: Should this be done here or before factorizing by checking that isMatrixModified? 
            var Kff = new Dictionary<int, IMatrixView>();
            foreach (int s in linearSystems.Keys) Kff[s] = linearSystems[s].Matrix;
            stiffnessDistribution.Update(Kff);
            subdomainGlobalMapping = new FetiDPSubdomainGlobalMapping(model, dofSeparator, stiffnessDistribution);
        }

        public Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix)
        {
            throw new NotImplementedException();
        }

        public void OrderDofs(bool alsoOrderConstrainedDofs) //TODO: Only order dofs of subdomains that are modified
        {
            var watch = new Stopwatch();
            watch.Start();

            // Order dofs
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (!subdomain.ConnectivityModified) continue; //TODO: Not sure about this

                matrixManagers[subdomain.ID].HandleDofOrderingWillBeModified();
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
                if (alsoOrderConstrainedDofs) subdomain.ConstrainedDofOrdering = dofOrderer.OrderConstrainedDofs(subdomain);

                // The next must done by the analyzer, so that subdomain.Forces is retained when doing back to back analyses.
                //subdomain.Forces = linearSystem.CreateZeroVector();
            }

            // Log dof statistics
            watch.Stop();
            Logger.LogTaskDuration("Dof ordering", watch.ElapsedMilliseconds);
            Logger.LogNumDofs("Global dofs", globalOrdering.NumGlobalFreeDofs);
        }

        public void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        public void Solve()
        {
            PrintLagrangeEqsData();
            PrintCornerEqsData();

            var watch = new Stopwatch();
            foreach (var linearSystem in linearSystems.Values)
            {
                if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            }

            // Separate the force vector
            watch.Start();
            var fr = new Dictionary<int, Vector>();
            var fbc = new Dictionary<int, Vector>();
            foreach (int s in subdomains.Keys)
            {
                int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                Vector f = linearSystems[s].RhsConcrete;
                fr[s] = f.GetSubvector(remainderDofs);
                fbc[s] = f.GetSubvector(cornerDofs);
            }
            watch.Stop();
            Logger.LogTaskDuration("Separating vectors & matrices", watch.ElapsedMilliseconds);
            watch.Reset();

            if (isStiffnessModified)
            {
                // Separate the stiffness matrix
                watch.Start();
                foreach (int s in subdomains.Keys)
                {
                    if (!subdomains[s].StiffnessModified) continue;
                    IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];
                    int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                    int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                    matrices.ExtractKrr(remainderDofs);
                    matrices.ExtractKcrKrc(cornerDofs, remainderDofs);
                    matrices.ExtractKcc(cornerDofs);
                }
                watch.Stop();
                Logger.LogTaskDuration("Separating vectors & matrices", watch.ElapsedMilliseconds);


                // Reorder internal dofs if needed by the preconditioner. TODO: Should I have done this previously in Initialize()?
                watch.Start();
                if (preconditionerFactory.ReorderInternalDofsForFactorization)
                {
                    foreach (ISubdomain subdomain in model.Subdomains)
                    {
                        if (subdomain.ConnectivityModified)
                        {
                            Debug.WriteLine($"{this.GetType().Name}: Reordering internal dofs of subdomain {subdomain.ID}.");
                            matrixManagers[subdomain.ID].ReorderInternalDofs(dofSeparator, subdomain);
                        }
                    }
                }

                // Calculate the preconditioner before factorizing each subdomain's Kff 
                preconditioner = preconditionerFactory.CreatePreconditioner(model, stiffnessDistribution, dofSeparator,
                    lagrangeEnumerator, matrixManagersGeneral);
                watch.Stop();
                Logger.LogTaskDuration("Calculating preconditioner", watch.ElapsedMilliseconds);

                // Factorize each subdomain's Krr
                watch.Restart();
                foreach (int s in subdomains.Keys)
                {
                    if (!subdomains[s].StiffnessModified) continue;
                    //TODO: If I can reuse Krr, I can also reuse its factorization. Therefore this must be inPlace. In contrast, FETI-1 needs Kff intact for Stiffness distribution, in the current design).
                    Debug.WriteLine(
                        $"{this.GetType().Name}: Inverting the remainder-remainder stiffness matrix of subdomain {s} in place.");
                    matrixManagers[s].InvertKrr(true);
                }
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);

                // Define FETI-DP flexibility matrices
                watch.Restart();
                flexibility = new FetiDPFlexibilityMatrix(dofSeparator, lagrangeEnumerator, matrixManagers);

#if DEBUG
                Matrix FIrr = MultiplyWithIdentity(lagrangeEnumerator.NumLagrangeMultipliers, lagrangeEnumerator.NumLagrangeMultipliers, flexibility.MultiplyFIrr);

#endif

                // Static condensation of remainder dofs (Schur complement).
                coarseProblemSolver.CreateAndInvertCoarseProblemMatrix(CornerNodesOfSubdomains, dofSeparator, matrixManagers);
                watch.Stop();
                Logger.LogTaskDuration("Setting up interface problem", watch.ElapsedMilliseconds);
                watch.Reset();

                isStiffnessModified = false;
            }

            // Static condensation for the force vectors
            watch.Start();
            Vector globalFcStar = coarseProblemSolver.CreateCoarseProblemRhs(dofSeparator, matrixManagers, fr, fbc);

            // Calculate the rhs vectors of the interface system
            Vector dr = CalcDisconnectedDisplacements(fr);
            double globalForcesNorm = CalcGlobalForcesNorm();
            watch.Stop();
            Logger.LogTaskDuration("Setting up interface problem", watch.ElapsedMilliseconds);

            // Solve the interface problem
            watch.Restart();
            (Vector lagranges, Vector uc) = interfaceProblemSolver.SolveInterfaceProblem(flexibility, preconditioner, 
                coarseProblemSolver, globalFcStar, dr, globalForcesNorm, Logger);
            watch.Stop();
            Logger.LogTaskDuration("Solving interface problem", watch.ElapsedMilliseconds);

            // Calculate the displacements of each subdomain
            watch.Restart();
            Dictionary<int, Vector> actualDisplacements = CalcActualDisplacements(lagranges, uc,  fr);
            foreach (var idSystem in linearSystems) idSystem.Value.SolutionConcrete = actualDisplacements[idSystem.Key];
            watch.Stop();
            Logger.LogTaskDuration("Calculate displacements from lagrange multipliers", watch.ElapsedMilliseconds);

            Logger.IncrementAnalysisStep();
        }

        private void PrintLagrangeEqsData()
        {
            List<int>[] subdomainLocals = new List<int>[subdomains.Count()];
            List<int>[] subdomainGobEqs = new List<int>[subdomains.Count()];
            List<int>[] subdomainPlMinus = new List<int>[subdomains.Count()];

            for (int i1 = 0; i1 < model.Subdomains.Count(); i1++)
            {
                subdomainLocals[i1] = new List<int>();
                subdomainGobEqs[i1] = new List<int>();
                subdomainPlMinus[i1] = new List<int>();
            }

            for (int i = 0; i < lagrangeEnumerator.NumLagrangeMultipliers; i++)
            {
                var lagrange = lagrangeEnumerator.LagrangeMultipliers[i];
                int subdId_pl = lagrange.SubdomainPlus.ID;
                int subdId_mn = lagrange.SubdomainMinus.ID;
                int dofId_local_plus = lagrange.SubdomainPlus.FreeDofOrdering.FreeDofs[lagrange.Node, lagrange.DofType];
                int dofId_local_minus = lagrange.SubdomainMinus.FreeDofOrdering.FreeDofs[lagrange.Node, lagrange.DofType];

                subdomainLocals[subdId_pl].Add(dofId_local_plus+1);
                subdomainGobEqs[subdId_pl].Add(i+1);
                subdomainPlMinus[subdId_pl].Add(+1);

                subdomainLocals[subdId_mn].Add(dofId_local_minus + 1);
                subdomainGobEqs[subdId_mn].Add(i + 1);
                subdomainPlMinus[subdId_mn].Add(-1);
            }

            var LagrangeNodeIds = new List<int>();
            var LagrangeNodeCouplingWaysNums = new Dictionary<int,int>();
            int previousNodeId = -1;
            int position = 0;
            for (int i = 0; i < lagrangeEnumerator.NumLagrangeMultipliers; i++)
            {
                int currentNodeId = lagrangeEnumerator.LagrangeMultipliers[i].Node.ID;
                if (!(currentNodeId==previousNodeId))
                {
                    LagrangeNodeIds.Add(currentNodeId);
                    position = LagrangeNodeIds.Count();
                    LagrangeNodeCouplingWaysNums[position] = 1;
                    previousNodeId = currentNodeId;
                    

                }
                else
                {
                    LagrangeNodeCouplingWaysNums[position] = LagrangeNodeCouplingWaysNums[position]+1;
                    //coupling_ways_num += 1;
                }
            }



            //string subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug\RVE_database\rve_no_1";
            string subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug_corner\RVE_database\rve_no_1";


            for (int i1 = 0; i1 < model.Subdomains.Count(); i1++)
            {
                var path1 = subdomainOutputPath_gen + @"\subdomainLagrangesLocals" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainLocals[i1], path1);

                path1 = subdomainOutputPath_gen + @"\subdomainLagrangesGlobalEqs" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainGobEqs[i1], path1);

                path1 = subdomainOutputPath_gen + @"\subdomainLagrangesPlMinus" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainPlMinus[i1], path1);

                path1 = subdomainOutputPath_gen + @"\RB_Nodes_IDs_MSOLVE_wise"+".txt";
                WriteToFileVector(LagrangeNodeIds, path1);

                path1 = subdomainOutputPath_gen + @"\RB_Nodes_IDs_couplingwaysNum_MSOLVE_wise"+".txt";
                WriteToFileVector(LagrangeNodeCouplingWaysNums.Values.ToList(), path1);
            }


        }

        private void PrintCornerEqsData()
        {
            List<int>[] subdomainLocals = new List<int>[subdomains.Count()];
            List<int>[] subdomainGobEqs = new List<int>[subdomains.Count()];
            
            for (int i1 = 0; i1 < model.Subdomains.Count(); i1++)
            {
                subdomainLocals[i1] = new List<int>();
                subdomainGobEqs[i1] = new List<int>();

                var subdomain = model.Subdomains[i1];

                foreach ((INode node, IDofType dof, _ ) in dofSeparator.SubdomainCornerDofOrderings[subdomain.ID])
                {
                    subdomainLocals[i1].Add(subdomain.FreeDofOrdering.FreeDofs[node, dof]+1);
                    subdomainGobEqs[i1].Add(dofSeparator.GlobalCornerDofOrdering[node, dof]+1);
                }
            }

            //string subdomainOutputPath_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug\RVE_database\rve_no_1";
            string subdomainOutputPath_gen = @"C:\Users\acivi\Documents\notes_elegxoi\REFERENCE_kanonikh_gewmetria_fe2_post_dg\REF2_10__000_renu_new_multiple_algorithms_check_develop_gia_fe2_3grsh_4182dofs_multiple2b_debug_corner\RVE_database\rve_no_1";


            for (int i1 = 0; i1 < model.Subdomains.Count(); i1++)
            {
                var path1 = subdomainOutputPath_gen + @"\subdomainCornerLocals" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainLocals[i1], path1);

                path1 = subdomainOutputPath_gen + @"\subdomainCornerGlobalEqs" + (i1 + 1) + ".txt";
                WriteToFileVector(subdomainGobEqs[i1], path1);
                
            }

        }

        public static void WriteToFileVector(List<int> array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.Count(); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();

            writer2.Dispose();

        }

        private static Matrix MultiplyWithIdentity(int numRows, int numCols, Action<Vector, Vector> matrixVectorMultiplication)
        {
            var result = Matrix.CreateZero(numRows, numCols);
            for (int j = 0; j < numCols; ++j)
            {
                var lhs = Vector.CreateZero(numCols);
                lhs[j] = 1.0;
                var rhs = Vector.CreateZero(numRows);
                matrixVectorMultiplication(lhs, rhs);
                result.SetSubcolumn(j, rhs);
            }
            return result;
        }

        /// <summary>
        /// Does not mutate this object.
        /// </summary>
        internal Dictionary<int, Vector> CalcActualDisplacements(Vector lagranges, Vector cornerDisplacements, 
            Dictionary<int, Vector> fr)
        {
            var freeDisplacements = new Dictionary<int, Vector>();
            foreach (int s in subdomains.Keys)
            {
                IFetiDPSubdomainMatrixManager matrices = matrixManagers[s];

                // ur[s] = inv(Krr[s]) * (fr[s] - Br[s]^T * lagranges - Krc[s] * Lc[s] * uc)
                Vector BrLambda = lagrangeEnumerator.BooleanMatrices[s].Multiply(lagranges, true);
                Vector KrcLcUc = dofSeparator.CornerBooleanMatrices[s].Multiply(cornerDisplacements);
                KrcLcUc = matrices.MultiplyKrcTimes(KrcLcUc);
                Vector temp = fr[s].Copy();
                temp.SubtractIntoThis(BrLambda);
                temp.SubtractIntoThis(KrcLcUc);
                Vector ur = matrices.MultiplyInverseKrrTimes(temp);

                // uf[s] = union(ur[s], ubc[s])
                // Remainder dofs
                var uf = Vector.CreateZero(subdomains[s].FreeDofOrdering.NumFreeDofs);
                int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                uf.CopyNonContiguouslyFrom(remainderDofs, ur);

                // Corner dofs: ubc[s] = Bc[s] * uc
                Vector ubc = dofSeparator.CornerBooleanMatrices[s].Multiply(cornerDisplacements);
                int[] cornerDofs = dofSeparator.CornerDofIndices[s];
                uf.CopyNonContiguouslyFrom(cornerDofs, ubc);

                freeDisplacements[s] = uf;
            }
            return freeDisplacements;
        }

        /// <summary>
        /// d = sum(Bs * generalInverse(Ks) * fs), where fs are the nodal forces applied to the dofs of subdomain s.
        /// Does not mutate this object.
        /// </summary>
        internal Vector CalcDisconnectedDisplacements(Dictionary<int, Vector> fr)
        {
            // dr = sum_over_s( Br[s] * inv(Krr[s]) * fr[s])
            var dr = Vector.CreateZero(lagrangeEnumerator.NumLagrangeMultipliers);
            foreach (int s in linearSystems.Keys)
            {
                SignedBooleanMatrixColMajor Br = lagrangeEnumerator.BooleanMatrices[s];
                Vector temp = matrixManagers[s].MultiplyInverseKrrTimes(fr[s]);
                temp = Br.Multiply(temp);
                dr.AddIntoThis(temp);
            }
            return dr;
        }

        /// <summary>
        /// Calculate the norm of the forces vector |f| = |K*u|. It is needed to check the convergence of PCG/PCPG.
        /// </summary>
        private double CalcGlobalForcesNorm()
        {
            //TODO: It would be better to do that using the global vector to avoid the homogeneous/heterogeneous averaging
            //      That would require the analyzer to build the global vector too though. Caution: we cannot take just 
            //      the nodal loads from the model, since the effect of other loads is only taken into account int 
            //      linearSystem.Rhs. Even if we could easily create the global forces vector, it might be wrong since 
            //      the analyzer may modify some of these loads, depending on time step, loading step, etc.
            var subdomainForces = new Dictionary<int, IVectorView>();
            foreach (var linearSystem in linearSystems.Values)
            {
                subdomainForces[linearSystem.Subdomain.ID] = linearSystem.RhsConcrete;
            }
            return subdomainGlobalMapping.CalculateGlobalForcesNorm(subdomainForces);
        }

        public class Builder
        {
            private readonly ICornerNodeSelection cornerNodeSelection;
            private readonly IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory; 

            public Builder(ICornerNodeSelection cornerNodeSelection, 
                IFetiDPSubdomainMatrixManagerFactory matrixManagerFactory)
            {
                this.cornerNodeSelection = cornerNodeSelection;
                this.matrixManagerFactory = matrixManagerFactory;
            }

            //TODO: We need to specify the ordering for remainder and possibly internal dofs, while IDofOrderer only works for free dofs.
            public IDofOrderer DofOrderer { get; set; } =
                new ReusingDofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public IFetiDPInterfaceProblemSolver InterfaceProblemSolver { get; set; } 
                = new FetiDPInterfaceProblemSolver.Builder().Build();
            public IFetiPreconditionerFactory PreconditionerFactory { get; set; } = new LumpedPreconditioner.Factory();
            public bool ProblemIsHomogeneous { get; set; } = true;

            public FetiDPSolverPrint BuildSolver(IStructuralModel model)
                => new FetiDPSolverPrint(model, cornerNodeSelection, matrixManagerFactory, DofOrderer, PreconditionerFactory,
                        ProblemIsHomogeneous, InterfaceProblemSolver);            
        }
    }
}
