using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Solvers.Direct;
using System;
using System.IO;
using System.Linq;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.Augmentation;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.StiffnessMatrices;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d.StiffnessMatrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP3d;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Analyzers.Loading;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Problems;

namespace ISAAR.MSolve.Tests
{
    public static class ShellPlateFetiDPExample
    {

        static double load = 0.0001;
        static bool paktwsi = true;
        [Fact]
        public static void RunExample()
        {
            CnstValues.preventOutputFileWrite();
            int subdiscrShell = 4; int discrShell = 3;
            (Model model, UsedDefinedCornerNodes cornerNodeSelection, UserDefinedMidsideNodes midSideNodeSelection, int loadedNodeId) = CreateModel(subdiscrShell, discrShell);
            double solutionz_Feti = Solve(model, cornerNodeSelection, midSideNodeSelection, loadedNodeId) ;

            (Model modelSerial, UsedDefinedCornerNodes corners, UserDefinedMidsideNodes midside, int loadedNodeIdserial) = CreateModelSerial(subdiscrShell, discrShell);
            double solutionz_serial = SolveSerial(modelSerial, corners, midside, loadedNodeIdserial);


            double sol2 = (solutionz_serial-solutionz_Feti)/solutionz_serial;
        }

        private static double Solve(Model model, UsedDefinedCornerNodes cornerNodeSelection, UserDefinedMidsideNodes midSideNodeSelection, int loadedNodeId)
        {
            var pcgSettings = new PcgSettings()
            {
                ConvergenceTolerance = 1E-4,
                MaxIterationsProvider = new FixedMaxIterationsProvider(1000)
            };
            var fetiMatrices = new FetiDPMatrixManagerFactorySkyline(new OrderingAmdSuiteSparse());
            //var fetiMatrices = new SkylineFetiDPSubdomainMatrixManager.Factory();
            //var fetiMatrices = new DenseFetiDPSubdomainMatrixManager.Factory();

            var fetiSolverBuilder = new FetiDPSolverSerial.Builder(fetiMatrices);  //A.3
            //var matrixManagerFactory = new FetiDP3dMatrixManagerFactoryDense();   //A.3.1
            //var matrixManagerFactory = new FetiDP3dMatrixManagerFactorySkyline();   //A.3.1
            //var matrixManagerFactory = new FetiDP3dMatrixManagerFactorySuiteSparse(new OrderingAmdSuiteSparse());   //A.3.1
            //var fetiSolverBuilder = new FetiDP3dSolverSerial.Builder(matrixManagerFactory);  //A.3

            //fetiSolverBuilder.InterfaceProblemSolver = interfaceSolverBuilder.Build();
            fetiSolverBuilder.StiffnessDistribution = StiffnessDistributionType.HeterogeneousCondensed;
            fetiSolverBuilder.Preconditioning = new DirichletPreconditioning();
            fetiSolverBuilder.PcgSettings = pcgSettings;

            ICrosspointStrategy crosspointStrategy = new MinimumConstraints(); //A.2

            fetiSolverBuilder.CrosspointStrategy = crosspointStrategy;

            FetiDPSolverSerial fetiSolver = fetiSolverBuilder.Build(model, cornerNodeSelection); //A.1
            //FetiDP3dSolverSerial fetiSolver = fetiSolverBuilder.Build(model, cornerNodeSelection, midSideNodeSelection); //A.1

            RunAnalysis(model, fetiSolver);

            var solution = fetiSolver.GatherGlobalDisplacements();

            int loadedDofzId = model.GlobalDofOrdering.GlobalFreeDofs[model.NodesDictionary[loadedNodeId], StructuralDof.TranslationZ];
            double solutionDofZ = solution[loadedDofzId];

            int subdID = model.NodesDictionary[loadedNodeId].SubdomainsDictionary.ElementAt(0).Value.ID;
            int localSubdomainDofId = model.SubdomainsDictionary[subdID].FreeDofOrdering.FreeDofs[model.NodesDictionary[loadedNodeId], StructuralDof.TranslationZ];
            double subdlocalSolutionZ = fetiSolver.GetLinearSystem(model.SubdomainsDictionary[subdID]).Solution[localSubdomainDofId];

            if(Math.Abs(solutionDofZ-subdlocalSolutionZ)/Math.Abs(solutionDofZ)>1e-6)
            {
                throw new Exception();
            }

            return subdlocalSolutionZ;


        }

        private static double SolveSerial(Model model, UsedDefinedCornerNodes cornerNodeSelection, UserDefinedMidsideNodes midSideNodeSelection, int loadedNodeId)
        {
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            // Analyzers
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
            //NewmarkDynamicAnalyzer parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);

            // Request output
            //int[] outputPositions = new int[model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs];
            //for( int i1 = 0;  i1 < model.SubdomainsDictionary[0].FreeDofOrdering.NumFreeDofs ; i1++ )
            //{

            //}
            //childAnalyzer.LogFactories[0] = new LinearAnalyzerLogFactory(new int[] { 0 });

            // Run the anlaysis - part a
            parentAnalyzer.Initialize();

           

            // Run the anlaysis - part b
            parentAnalyzer.Solve();


            var globalU = solver.LinearSystems[1].Solution.CopyToArray();//

            

            int loadedDofzId = model.GlobalDofOrdering.GlobalFreeDofs[model.NodesDictionary[loadedNodeId], StructuralDof.TranslationZ];
            double solutionDofZ = globalU[loadedDofzId];

            
            return solutionDofZ;

        }

        internal static void RunAnalysis(IModel model, ISolverMpi solver)
        {
            // Run the analysis
            solver.OrderDofs(false);
            foreach (ISubdomain subdomain in model.EnumerateSubdomains())
            {
                ILinearSystemMpi linearSystem = solver.GetLinearSystem(subdomain);
                linearSystem.Reset(); // Necessary to define the linear system's size 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }
            solver.BuildGlobalMatrix(new ElementStructuralStiffnessProvider());
            model.ApplyLoads();
            LoadingUtilities.ApplyNodalLoads(model, solver);
            solver.Solve();
        }

        public static (Model model, UsedDefinedCornerNodes cornerNodeSelection, UserDefinedMidsideNodes midSideNodeSelection, int loadedNodeId) 
            CreateModel(int subdiscrShell, int discrShell)
        {
             
            int elem1 = subdiscrShell * discrShell;
            int elem2 = subdiscrShell * discrShell;
            var mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(2, 2, 2, subdiscrShell, discrShell);
            var gp = mpgp.Item2;
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;

            var coefficientsProvider = new SpectralRepresentation2DRandomField(0, 0, 0, 0.01, 2 * Math.PI / ((double)20), 20);
            coefficientsProvider.RandomVariables = null;
            double[] o_xsunol = FEMMeshBuilder.ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, gp.L1, gp.L2, gp.elem1, gp.elem2, gp.a1_shell, new double[] { 0, 0, 0 },
               new o_x_parameters() , coefficientsProvider);


            //model and node coordinates
            Model model = new Model();
            int totalNodes = o_xsunol.Length / 6;
            int[] newNumbering = new int[totalNodes];
            for(int i1 = 0; i1 < totalNodes; i1++) { newNumbering[i1] = i1 + 1; }
           
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int eswterikosNodeCounter = 0;
            var renumbering = new renumbering(newNumbering);
            for(int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter +  1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;

            //material
            var material2 = new ShellElasticMaterial3D()
            {
                YoungModulus = gp.E_shell,
                PoissonRatio = gp.ni_shell,
                ShearCorrectionCoefficientK = 5 / 6,
            };


            //perioxh orismou twn elements
            #region
            int elements = elem1 * elem2;
            int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
            int komvoi_8 = fdof_8 / 5;
            int[,] t_shell;
            t_shell = topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

            double[] Tk_vec = new double[8];
            double[][] VH = new double[8][];
            int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
            Element e2;
            int ElementID;

            int[] subdIDs = new int[discrShell * discrShell];
            for (int i1 = 0; i1 < discrShell*discrShell; i1++)
            {
                model.SubdomainsDictionary.Add(i1+1, new Subdomain(i1+1));
            }
            

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = gp.tk;
            }

            int eswterikosElementCounter = 0;
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter +  1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    //
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
                    {
                        //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] ), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] )]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                int subdomainID = GetNodeSubdomainID(discrShell, discrShell, e2, gp.L1, gp.L2);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            #endregion
            int loadedNodeId = AddConstraintsAndLoads(model,gp.L1, gp.L2);

            model.ConnectDataStructures();

            List<Node> cornerNodeList = model.NodesDictionary.Values.Where(x => x.SubdomainsDictionary.Keys.Count == 4).ToList();

            var cornerNodes = DefineCornerNodesPerSubdomainAndOtherwise(cornerNodeList, model);

            var cornerNodes_ = cornerNodes.Select(x => ((ISubdomain)model.SubdomainsDictionary[x.Key], x.Value)).ToDictionary(x => x.Item1, x => x.Value);

            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes_);

            List<List<int>> ExtraconstraintNodes = GetExtraConstraintsOfplate(discrShell, subdiscrShell, gp.L1, gp.L2, model);

            Dictionary<ISubdomain, HashSet<INode>> extraConstrNodesofsubd = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subd in model.EnumerateSubdomains()) extraConstrNodesofsubd.Add((ISubdomain)subd, new HashSet<INode>());
            foreach (var extraConstrNodeList in ExtraconstraintNodes)
            {
                int extraNodeId = extraConstrNodeList[0];
                foreach (var subd in model.NodesDictionary[extraNodeId].SubdomainsDictionary.Values)
                {
                    extraConstrNodesofsubd[subd].Add(model.NodesDictionary[extraNodeId]);
                }
            }

            var midSideNodeSelection = new UserDefinedMidsideNodes(extraConstrNodesofsubd,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });

            
            return (model, cornerNodeSelection, midSideNodeSelection, loadedNodeId);
        }

        public static (Model model, UsedDefinedCornerNodes cornerNodeSelection, UserDefinedMidsideNodes midSideNodeSelection, int loadedNodeId)
            CreateModelSerial(int subdiscrShell, int discrShell)
        {
            int subdId = 1;
            int elem1 = subdiscrShell * discrShell;
            int elem2 = subdiscrShell * discrShell;
            var mpgp = GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(2, 2, 2, subdiscrShell, discrShell);
            var gp = mpgp.Item2;
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;

            var coefficientsProvider = new SpectralRepresentation2DRandomField(0, 0, 0, 0.01, 2 * Math.PI / ((double)20), 20);
            coefficientsProvider.RandomVariables = null;
            double[] o_xsunol = FEMMeshBuilder.ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, gp.L1, gp.L2, gp.elem1, gp.elem2, gp.a1_shell, new double[] { 0, 0, 0 },
               new o_x_parameters(), coefficientsProvider);


            //model and node coordinates
            Model model = new Model();
            int totalNodes = o_xsunol.Length / 6;
            int[] newNumbering = new int[totalNodes];
            for (int i1 = 0; i1 < totalNodes; i1++) { newNumbering[i1] = i1 + 1; }

            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int eswterikosNodeCounter = 0;
            var renumbering = new renumbering(newNumbering);
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;

            //material
            var material2 = new ShellElasticMaterial3D()
            {
                YoungModulus = gp.E_shell,
                PoissonRatio = gp.ni_shell,
                ShearCorrectionCoefficientK = 5 / 6,
            };


            //perioxh orismou twn elements
            #region
            int elements = elem1 * elem2;
            int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
            int komvoi_8 = fdof_8 / 5;
            int[,] t_shell;
            t_shell = topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

            double[] Tk_vec = new double[8];
            double[][] VH = new double[8][];
            int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
            Element e2;
            int ElementID;

            
             model.SubdomainsDictionary.Add(1, new Subdomain(1));
            


            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = gp.tk;
            }

            int eswterikosElementCounter = 0;
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    //
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
                    {
                        //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1]), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1])]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[1].Elements.Add(e2.ID, e2);
                eswterikosElementCounter++;
            }
            #endregion
            int loadedNodeId = AddConstraintsAndLoads(model, gp.L1, gp.L2);

            model.ConnectDataStructures();

            List<Node> cornerNodeList = model.NodesDictionary.Values.Where(x => x.SubdomainsDictionary.Keys.Count == 4).ToList();

            var cornerNodes = DefineCornerNodesPerSubdomainAndOtherwise(cornerNodeList, model);

            var cornerNodes_ = cornerNodes.Select(x => ((ISubdomain)model.SubdomainsDictionary[x.Key], x.Value)).ToDictionary(x => x.Item1, x => x.Value);

            var cornerNodeSelection = new UsedDefinedCornerNodes(cornerNodes_);

            List<List<int>> ExtraconstraintNodes = GetExtraConstraintsOfplate(discrShell, subdiscrShell, gp.L1, gp.L2, model);

            Dictionary<ISubdomain, HashSet<INode>> extraConstrNodesofsubd = new Dictionary<ISubdomain, HashSet<INode>>();
            foreach (Subdomain subd in model.EnumerateSubdomains()) extraConstrNodesofsubd.Add((ISubdomain)subd, new HashSet<INode>());
            foreach (var extraConstrNodeList in ExtraconstraintNodes)
            {
                int extraNodeId = extraConstrNodeList[0];
                foreach (var subd in model.NodesDictionary[extraNodeId].SubdomainsDictionary.Values)
                {
                    extraConstrNodesofsubd[subd].Add(model.NodesDictionary[extraNodeId]);
                }
            }

            var midSideNodeSelection = new UserDefinedMidsideNodes(extraConstrNodesofsubd,
                new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });


            return (model, cornerNodeSelection, midSideNodeSelection, loadedNodeId);
        }

        private static int AddConstraintsAndLoads(Model model, double L1, double L2)
        {
            int loadedNodeId=-1; 
            double tol = 0.00001;
            foreach (Node node in model.NodesDictionary.Values)
            {
                if(Math.Abs(node.X-0.5*L1)<tol| Math.Abs(node.X + 0.5 * L1) < tol | Math.Abs(node.Y + 0.5 * L1) < tol | Math.Abs(node.Y - 0.5 * L1) < tol)
                {
                    node.Constraints.Add(new Discretization.Constraint() { DOF = StructuralDof.TranslationX, Amount = 0 });
                    node.Constraints.Add(new Discretization.Constraint() { DOF = StructuralDof.TranslationY, Amount = 0 });
                    node.Constraints.Add(new Discretization.Constraint() { DOF = StructuralDof.TranslationZ, Amount = 0 });
                    if (paktwsi)
                    {
                        node.Constraints.Add(new Discretization.Constraint() { DOF = StructuralDof.RotationX, Amount = 0 });
                        node.Constraints.Add(new Discretization.Constraint() { DOF = StructuralDof.RotationY, Amount = 0 });
                    }
                }
                if (Math.Abs(node.X - 0) < tol &&  Math.Abs(node.Y -0) < tol )
                {
                    model.Loads.Add(new Load() { Amount = load, Node = node, DOF = StructuralDof.TranslationZ });
                    loadedNodeId = node.ID;
                }

            }

            return loadedNodeId;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 95, //150, // diastaseis
                L02 = 95, //150,
                L03 = 95, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr1 * subdiscr1,
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = discr1_shell * subdiscr1_shell,
                elem2 = discr1_shell * subdiscr1_shell,
                L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 50,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.20, //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.25, //0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.20, //0.05,// Gpa
                D_o_1 = 0.25, //0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static int[,] topologia_shell_coh(int elements, int elem1, int elem2, object komvoi_8)
        {
            int elem;
            int[,] t_shell = new int[elements, 8];
            for (int nrow = 0; nrow < elem1; nrow++)
            {
                for (int nline = 0; nline < elem2; nline++)
                {
                    elem = (nrow + 1 - 1) * elem2 + nline + 1;//nrow+ 1 nline+1 einai zero based 
                    t_shell[elem - 1, -1 + 1] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 8] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 4] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;

                    t_shell[elem - 1, -1 + 5] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 2;
                    t_shell[elem - 1, -1 + 7] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 1;

                    t_shell[elem - 1, -1 + 2] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 6] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 3] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;
                }
            }
            return t_shell;

        }

        public static int GetNodeSubdomainID(int n1, int n2, Element element, double L01, double L02)
        {
            var x_meso = element.NodesDictionary.Values.Select(x => x.X).ToArray().Average() +0.5*L01;
            var y_meso = element.NodesDictionary.Values.Select(x => x.Y).ToArray().Average() +0.5*L02;

            int x_gemata = (int)Math.Truncate(x_meso/(L01 / n1));
            int y_gemata = (int)Math.Truncate(y_meso / (L02 / n2));

            int subd_ID = n2 * x_gemata + y_gemata + 1;
            return subd_ID;

        }

        private  static Dictionary<int, HashSet<INode>> DefineCornerNodesPerSubdomainAndOtherwise(List<Node>  CornerNodesList, Model model)
        {
            // a copy of this is used in modelCreatorInput class
            Dictionary<int, HashSet<INode>> cornerNodesList = new Dictionary<int, HashSet<INode>>(model.EnumerateSubdomains().Count());
            
            foreach (Subdomain subdomain in model.EnumerateSubdomains())
            {
                cornerNodesList.Add(subdomain.ID, new HashSet<INode>());
            }

            foreach (Node node1 in CornerNodesList)
            {
                foreach (Subdomain subdomain in node1.SubdomainsDictionary.Values)
                {
                    cornerNodesList[subdomain.ID].Add(node1);
                }
            }

            //foreach (Subdomain subdomain in model.Subdomains)
            //{
            //    cornerNodes.Add(subdomain.ID, cornerNodesList[subdomain.ID]);
            //}

            return cornerNodesList;
        }

        private static List<List<int>> GetExtraConstraintsOfplate(int discrShell, int subdiscrShell, double l1, double l2, Model model)
        {
            if (!(l1 == l2)) { throw new NotImplementedException(); }
            List<List<int>> extraconstraintNodes = new List<List<int>>();

            double elem_length = l1 / (discrShell * subdiscrShell);
            double subd_length = l1 / discrShell;
            double subd_half_length = 0.5 * subd_length;
            double tol = 1e-8;
            double elem_half_length =0.5* l1 / (discrShell * subdiscrShell);


            foreach (Node node in model.NodesDictionary.Values)
            {
                if(node.SubdomainsDictionary.Values.Count==2)
                {
                    var x_coord = node.X + 0.5 * l1 + 0.0001 * elem_length;
                    bool isNotSubdBound = true;
                    bool isSubdomainHalf_x = false;

                    double x_round_subd = Math.Truncate(x_coord / subd_length) * subd_length;
                    double x_round_half_sub = Math.Truncate(x_coord / subd_half_length) * subd_half_length;
                    double x_elem_round = Math.Truncate(x_coord / elem_half_length) * elem_half_length;

                    if(Math.Abs(x_round_subd - x_elem_round) <tol) { isNotSubdBound = false; }
                    if (((x_round_half_sub - x_elem_round) < tol)&& isNotSubdBound) 
                    {
                        isSubdomainHalf_x = true; 
                        List<int> extraNodes = new List<int>() { node.ID };
                        extraconstraintNodes.Add(extraNodes);
                        continue;
                    }

                    var y_coord = node.Y + 0.5 * l1 + 0.0001 * elem_length;
                    bool isNotSubdBoundy = true;

                    double y_round_subd = Math.Truncate(y_coord / subd_length) * subd_length;
                    double y_round_half_sub = Math.Truncate(y_coord / subd_half_length) * subd_half_length;
                    double y_elem_round = Math.Truncate(y_coord / elem_half_length) * elem_half_length;

                    if (Math.Abs(y_round_subd - y_elem_round) < tol) { isNotSubdBoundy = false; }
                    if (((y_round_half_sub - y_elem_round) < tol) && isNotSubdBoundy)
                    {
                        List<int> extraNodes = new List<int>() { node.ID };
                        extraconstraintNodes.Add(extraNodes);
                        continue;
                    }
                }
            }

            return extraconstraintNodes;
        }

    }
}
