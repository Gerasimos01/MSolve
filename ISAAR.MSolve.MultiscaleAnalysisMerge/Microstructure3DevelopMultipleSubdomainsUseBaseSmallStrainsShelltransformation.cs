﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;

using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformation : StructuralProblemsMicrostructureBase, IShellMaterial //A.1
    {
        // proelefsi: Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrains2D
        // allages: transformation apo shell2D klp wste na kanei implement to Ishell
        
        //allagesII:i) comment out duplicate calculations of f2_vectors
        //ii) tlk to transformation eixe hdh perasthei

        //A.2
        public double[] NormalVectorV3 { get; set; }
        public double[] TangentVectorV1 { get; set; }
        public double[] TangentVectorV2 { get; set; }

        //TODO: create base class and use it for different microstructures-scale transitions.
        private Model model { get; set; }
        //private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private Dictionary<int, Node> boundaryNodes { get; set; }
        private IdegenerateRVEbuilder rveBuilder;
        //private NewtonRaphsonNonLinearAnalyzer microAnalyzer;
        private double volume;
        Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain;
        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements;
        private IScaleTransitions scaleTransitions = new SmallStrain3Dto2DplaneStressScaleTransition(); //TODO: mporoume na to dinoume ston constructor


        // aparaithta gia to implementation tou IFiniteElementMaterial3D
        Matrix2D constitutiveMatrix;
        private double[] trueStressVec; // TODO: rename stresses 
        Matrix2D transformationMatrix; // gia to shell
        private bool modified; // opws sto MohrCoulomb gia to modified

        private double[,] Cijrs_prev;
        private bool matrices_not_initialized = true;
        private double tol;
        public void InitializeMatrices()
        {
            Cijrs_prev = new double[3,3];
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
            constitutiveMatrix = new Matrix2D(new double[3,3]);
            this.CalculateTransformationMatrix(new Vector(TangentVectorV1), new Vector(TangentVectorV2));
        }
        
        //double[] Stresses { get; }
        //IMatrix2D ConstitutiveMatrix { get; } TODOGerasimos

        public Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformation(IdegenerateRVEbuilder rveBuilder)
        {
            this.rveBuilder = rveBuilder;
            Tuple<Model, Dictionary<int, Node>,double> modelAndBoundaryNodes = this.rveBuilder.GetModelAndBoundaryNodes();
            this.model = modelAndBoundaryNodes.Item1;
            this.boundaryNodes = modelAndBoundaryNodes.Item2;
            this.volume = modelAndBoundaryNodes.Item3;
            DefineAppropriateConstraintsForBoundaryNodes();
            this.model.ConnectDataStructures();
            this.InitializeFreeAndPrescribedDofsInitialDisplacementVectors();
        }

        private void DefineAppropriateConstraintsForBoundaryNodes()
        {
            foreach(Node boundaryNode in boundaryNodes.Values)
            {
                var RigidBodyNodeConstraints = rveBuilder.GetModelRigidBodyNodeConstraints(model);
                scaleTransitions.ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(model, boundaryNode, RigidBodyNodeConstraints);
            }
        }

        private void InitializeFreeAndPrescribedDofsInitialDisplacementVectors()
        {
            uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, Vector>();
            foreach(Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                uInitialFreeDOFDisplacementsPerSubdomain.Add(subdomain.ID, new Vector(subdomain.TotalDOFs));// prosoxh sto Id twn subdomain
            }            
            double[] smallStrainVec = new double[3] ;
            initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, initialConvergedBoundaryDisplacements);
            }            
        }

        public IShellMaterial Clone()
        {
            return new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformation(rveBuilder);
        }

        object ICloneable.Clone()
        {
            return new Microstructure3DevelopMultipleSubdomainsUseBaseSmallStrainsShelltransformation(rveBuilder);
        }

        public Dictionary<int, Node> BoundaryNodesDictionary
        {
            get { return boundaryNodes; }
        }
        public IList<Node> BoundaryNodes
        {
            get { return boundaryNodes.Values.ToList<Node>(); }
        }

        public void UpdateMaterial(double[] smallStrainVec)
        {

            if (matrices_not_initialized) 
            { this.InitializeMatrices(); } // mporei na mpei kai ston constructor to initialization an exoume kai ta surface basis vectors ekei

            double[] rveCoordinatesSmallStrainVec = TransformStrains(smallStrainVec);
            smallStrainVec = rveCoordinatesSmallStrainVec;

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    Cijrs_prev[i1, j1] = this.constitutiveMatrix[i1, j1];}
            }

            #region Rve prescribed Dofs total DIsplacement Dictionary Creation (nessesary for NRNLAnalyzer)
            // epivolh metakinhsewn ston analyzer pou molis dhmiourghsame --> tha ginetai sto dictionary
            Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                smallStrainVec, totalPrescribedBoundaryDisplacements);
            }
            #endregion


            

            //var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            //linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            // TODO alternatively:

            var linearSystems = CreateNecessaryLinearSystems(model);
            //var linearSystems = new Dictionary<int, ILinearSystem>();
            //foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)" tou opoiu ta stoixeia ginontai access kai me ID
            //{
            //    linearSystems.Add(subdomain.ID, new SkylineLinearSystem(subdomain.ID, subdomain.Forces));// prosoxh sto Id twn subdomain
            //}

            var solver = GetAppropriateSolver(linearSystems);
            //var solver = new SolverSkyline(linearSystems[1]); //TODO this depends on the number of model.SubdomainsDictionary.Values

            #region Creation of nessesary analyzers for NRNLAnalyzer and Creation of Microstructure analyzer (NRNLdevelop temporarilly) and solution ;
            int increments = 1; int MaxIterations = 100; int IterationsForMatrixRebuild = 1;
            (NewtonRaphsonNonLinearAnalyzerDevelop microAnalyzer, ProblemStructural provider, ElementStructuralStiffnessProvider elementProvider,
                SubdomainGlobalMapping[] subdomainMappers) = AnalyzeMicrostructure(model, linearSystems, solver, increments, MaxIterations, IterationsForMatrixRebuild,
                totalPrescribedBoundaryDisplacements, initialConvergedBoundaryDisplacements, boundaryNodes, uInitialFreeDOFDisplacementsPerSubdomain);
            //ProblemStructural provider = new ProblemStructural(model, linearSystems); //B.1 apo edw Methodos Analyze microstructure
            ////var linearSystemsArray = new[] { linearSystems[1] }; //those depend on the number of model.SubdomainsDictionary.Values as well
            ////var subdomainUpdaters = new[] { new NonLinearSubdomainUpdaterWithInitialConditions(model.Subdomains[0]) };
            ////var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            //int totalSubdomains = model.Subdomains.Count;
            //var linearSystemsArray = new ILinearSystem[totalSubdomains]; 
            //var subdomainUpdaters = new NonLinearSubdomainUpdaterWithInitialConditions[totalSubdomains];
            //var subdomainMappers = new SubdomainGlobalMapping[totalSubdomains];    int counter = 0;
            //foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)"
            //{
            //    linearSystemsArray[counter] = linearSystems[subdomain.ID];
            //    subdomainUpdaters[counter] = new NonLinearSubdomainUpdaterWithInitialConditions(subdomain);
            //    subdomainMappers[counter] = new SubdomainGlobalMapping(subdomain);
            //    counter++;
            //}


            //ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
            //Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            ////equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //{
            //    equivalentContributionsAssemblers.Add(subdomain.ID, new EquivalentContributionsAssebler(subdomain, elementProvider));
            //}

            //var increments = 1; // microAnalyzer onomazetai o childAnalyzer
            //NewtonRaphsonNonLinearAnalyzerDevelop microAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelop(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
            //    boundaryNodes, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, equivalentContributionsAssemblers);
            //microAnalyzer.SetMaxIterations = 100;
            //microAnalyzer.SetIterationsForMatrixRebuild = 1;


            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, microAnalyzer, linearSystems);
            //parentAnalyzer.BuildMatrices();
            //parentAnalyzer.Initialize();
            //parentAnalyzer.Solve(); //B.2 ews edw Methodos Analyze microstructure
            #endregion

            #region update of free converged displacements vectors
            uInitialFreeDOFDisplacementsPerSubdomain = microAnalyzer.GetConvergedSolutionVectorsOfFreeDofs();// ousiastika to u pou twra taftizetai me to uPlusuu
            #endregion


            #region INTEGRATION stresses 
            Dictionary<int, Vector> du = microAnalyzer.GetConvergedIncrementalSolutionVectorsOfFreeDofs();
            //A.6
            //double[] FppReactionVector = SubdomainCalculations.CalculateFppReactionsVector(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes,
            //    uInitialFreeDOFDisplacementsPerSubdomain[model.Subdomains[0].ID], du[model.Subdomains[0].ID], initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, increments, increments);
            Dictionary<int, double[]> FppReactionVectorSubdomains = SubdomainCalculationsMultiple.CalculateFppReactionsVectorSubdomains(model, elementProvider, scaleTransitions, boundaryNodes,
                uInitialFreeDOFDisplacementsPerSubdomain, du, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, increments, increments);
            double[] FppReactionVector= SubdomainCalculationsMultiple.CombineMultipleSubdomainsStressesIntegrationVectorsIntoTotal(FppReactionVectorSubdomains);



            double[] DqFpp = SubdomainCalculations.CalculateDqFpp(FppReactionVector, scaleTransitions, boundaryNodes);

            trueStressVec = new double [DqFpp.Length];
            for (int i1 = 0; i1 < DqFpp.Length; i1++)
            { trueStressVec[i1]=(1 / volume) * DqFpp[i1]; }

            //SPK_vec = new double[6] { SPK_mat[0,0], SPK_mat[1,1], SPK_mat[2,2], SPK_mat[0,1], SPK_mat[1,2], SPK_mat[0,2] };
            //TODOna elegxthei h parapanw anadiataxh kai o pollaplasiasmos
            #endregion

            #region INTEGRATION constitutive Matrix
            //TODOGERASIMOS edw thelei NEWtonRaphson.BuildMatrices kai sigoura to solver. Initialize kapou edw
            provider.Reset(); // prosoxh xwris to provider.reset to BuildMatrices den tha vrei null ta ks kai den tha ta xanahtisei
            microAnalyzer.BuildMatrices();

            //A.1 //TODO: comment out old commands
            //double[][] KfpDq = SubdomainCalculations.CalculateKfreeprescribedDqMultiplications(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);
            Dictionary<int, double[][]> KfpDqSubdomains = SubdomainCalculationsMultiple.CalculateKfreeprescribedDqMultiplicationsMultiple(model, elementProvider, scaleTransitions, boundaryNodes);


            // TODO: replace provider.Reset(); microAnalyzer.BuildMatrices(); 
            //Dictionary<int, double[][]> KfpDqSubdomains =...; Dictionary<int, double[][]> KppDqVectorsSubdomains =...;
            //with the following two commands
            //var boundaryElements = GetSubdomainsBoundaryFiniteElementsDictionaries(model, boundaryNodes);
            //(Dictionary<int, double[][]> KfpDqSubdomains, Dictionary<int, double[][]> KppDqVectorsSubdomains) =
            //    SubdomainCalculationsSimultaneous.UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultiple(model,
            //    elementProvider, scaleTransitions, boundaryNodes, boundaryElements,linearSystems);


            ////calculate matrices for debug
            //(var constrainedNodalDOFsDictionary,var  TotalConstrainedDOFs) = SubdomainCalculations.GetConstrainednodalDOFsDictionaryForDebugging(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);
            //var Dq = SubdomainCalculations.GetDqTransitionsMatrixForDebugging(model.Subdomains[0], elementProvider,  scaleTransitions,boundaryNodes,  constrainedNodalDOFsDictionary,  TotalConstrainedDOFs);
            //var Kfp = SubdomainCalculations.GetKfpMatrixForDebugging(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes, constrainedNodalDOFsDictionary, TotalConstrainedDOFs);
            //var Kff = SubdomainCalculations.CalculateGlobalMatrix(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes, constrainedNodalDOFsDictionary, TotalConstrainedDOFs);
            //var Kfp_Dq = (new Matrix2D(Kfp)) * (new Matrix2D(Dq));
            //var Kpp = SubdomainCalculations.GetKppMatrixForDebugging(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes, constrainedNodalDOFsDictionary, TotalConstrainedDOFs);
            //var KffinvKfpDq = (new Matrix2D(Kff)).SolveLU(Kfp_Dq, true);
            //var KpfKffinvKfpDq = (new Matrix2D(Kfp)).Transpose() * KffinvKfpDq;
            //var Kpp_Dq= (new Matrix2D(Kpp)) * (new Matrix2D(Dq));
            //var ektosDq = new double[Kpp_Dq.Rows, Kpp_Dq.Columns]; // - KpfKffinvKfpDq;
            //for (int i1 = 0; i1 < Kpp_Dq.Rows; i1++){for (int i2 = 0; i2 < Kpp_Dq.Columns; i2++){ ektosDq[i1, i2] = Kpp_Dq[i1, i2] - KpfKffinvKfpDq[i1, i2]; }};
            //var Dq_CondDq = (new Matrix2D(Dq)).Transpose() * (new Matrix2D(ektosDq));
            //var d2W_dfdf_ch = new Matrix2D(Dq_CondDq.Data); d2W_dfdf_ch.Scale((1 / volume));
            ////

            // to BUildmatirces to exoume krathsei panw exw apo th sunarthsh tou f2
            //A.2
            //double[][] f2_vectors = SubdomainCalculations.CalculateKffinverseKfpDq(KfpDq, model, elementProvider, scaleTransitions, boundaryNodes, solver, linearSystems);            
            Dictionary<int, double[][]> f2_vectorsSubdomains = SubdomainCalculationsMultiple.CalculateKffinverseKfpDqSubdomains(KfpDqSubdomains, model, elementProvider, scaleTransitions, boundaryNodes, solver, linearSystems, subdomainMappers);

            //A.3
            // f3_vectors = SubdomainCalculations.CalculateKpfKffinverseKfpDq(f2_vectors, model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);
            Dictionary<int, double[][]> f3_vectorsSubdomains = SubdomainCalculationsMultiple.CalculateKpfKffinverseKfpDqSubdomains(f2_vectorsSubdomains, model, elementProvider, scaleTransitions, boundaryNodes);

            //A.4
            // KppDqVectors = SubdomainCalculations.CalculateKppDqMultiplications(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);
            Dictionary<int, double[][]> KppDqVectorsSubdomains = SubdomainCalculationsMultiple.CalculateKppDqSubdomainsMultiplications(model, elementProvider, scaleTransitions, boundaryNodes);

            //A.5
            double[][] f3_vectors = SubdomainCalculationsMultiple.CombineMultipleSubdomainsIntegrationVectorsIntoTotal(f3_vectorsSubdomains,scaleTransitions);
            double[][] KppDqVectors = SubdomainCalculationsMultiple.CombineMultipleSubdomainsIntegrationVectorsIntoTotal(KppDqVectorsSubdomains,scaleTransitions);

            double[][] f4_vectors = SubdomainCalculations.SubtractConsecutiveVectors(KppDqVectors, f3_vectors);
            double[,] DqCondDq = SubdomainCalculations.CalculateDqCondDq(f4_vectors, scaleTransitions, boundaryNodes);

            double[,] constitutiveMat = new double[DqCondDq.GetLength(0), DqCondDq.GetLength(1)];
            for (int i1 = 0; i1 < DqCondDq.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < DqCondDq.GetLength(1); i2++)
                {
                    constitutiveMat[i1, i2] = (1 / volume) * DqCondDq[i1, i2];
                }
            }
            #endregion

            #region update of prescribed converged displacements vectors;
            initialConvergedBoundaryDisplacements = totalPrescribedBoundaryDisplacements;
            #endregion

            #region constitutive tensors transformation methods
            // transformation gia to shell 
            (var transformedTrueStressVec, var transformedConstitutiveMat) =StressesAndConstitutiveMatrixTransformation(trueStressVec, constitutiveMat);
            #endregion

            this.constitutiveMatrix = new Matrix2D(transformedConstitutiveMat);
            trueStressVec = transformedTrueStressVec;

            //PrintMethodsForDebug(KfpDq, f2_vectors, f3_vectors, KppDqVectors, f4_vectors, DqCondDq, d2W_dfdf, Cijrs);
            this.modified = CheckIfConstitutiveMatrixChanged(); 
        }

        
        private void CalculateTransformationMatrix(Vector surfaceBasisVector1, Vector surfaceBasisVector2)
        {
            var auxMatrix1 = new Matrix2D(2, 2);  //auxMatrix: covariant metric coefficients gab
            auxMatrix1[0, 0] = surfaceBasisVector1.DotProduct(surfaceBasisVector1);
            auxMatrix1[0, 1] = surfaceBasisVector1.DotProduct(surfaceBasisVector2);
            auxMatrix1[1, 0] = surfaceBasisVector2.DotProduct(surfaceBasisVector1);
            auxMatrix1[1, 1] = surfaceBasisVector2.DotProduct(surfaceBasisVector2);
            (Matrix2D inverse, double det) = auxMatrix1.Invert2x2AndDeterminant(); //inverse: contravariant metric coefficients g_ab (ekthetis ta a,b)

            //Contravariant base vectors
            double[][] G_i = new double[2][];
            for (int i1 = 0; i1 < 2; i1++)
            {
                G_i[i1] = new double[3];
                for (int i2 = 0; i2 < 3; i2++)
                {
                    G_i[i1][i2] = inverse[i1, 0] * surfaceBasisVector1[i2] + inverse[i1, 1] * surfaceBasisVector2[i2];
                }
            }

            //Normalised covariant base vectors
            double[][] Ei = new double[2][];// to trito den xreiazetai

            Ei[0] = new double[3]; Array.Copy(surfaceBasisVector1.Data, Ei[0], 3);
            double G1_norm = surfaceBasisVector1.Norm;
            for (int i1 = 0; i1 < 3; i1++) { Ei[0][i1] = Ei[0][i1] / G1_norm; }

            double G2_dot_E1 = 0;
            for (int i1 = 0; i1 < 3; i1++) { G2_dot_E1 += surfaceBasisVector2[i1] * Ei[0][i1]; }

            double[] projection = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { projection[i1] = G2_dot_E1 * Ei[0][i1]; }

            Ei[1] = new double[3];
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = surfaceBasisVector2[i1] - projection[i1]; }
            double norm1 = new Vector(Ei[1]).Norm;
            for (int i1 = 0; i1 < 3; i1++) { Ei[1][i1] = Ei[1][i1] / norm1; }

            double[,] EiDOTG_j = new double[2, 2];

            for (int i1 = 0; i1 < 2; i1++)
            {
                for (int i2 = 0; i2 < 2; i2++)
                {
                    EiDOTG_j[i1, i2] = new Vector(Ei[i1]).DotProduct(new Vector(G_i[i2]));
                }
            }

             transformationMatrix = new Matrix2D(new double[3, 3] { {EiDOTG_j[0,0]*EiDOTG_j[0,0],EiDOTG_j[0,1]*EiDOTG_j[0,1],EiDOTG_j[0,0]*EiDOTG_j[0,1]  },
                 {EiDOTG_j[1,0]*EiDOTG_j[1,0],EiDOTG_j[1,1]*EiDOTG_j[1,1],EiDOTG_j[1,0]*EiDOTG_j[1,1]  },
                {2*EiDOTG_j[1,0]*EiDOTG_j[0,0],2*EiDOTG_j[1,1]*EiDOTG_j[0,1],EiDOTG_j[1,0]*EiDOTG_j[0,1]+EiDOTG_j[1,1]*EiDOTG_j[0,0]   } });
        }

        private double[] TransformStrains(double[] smallStrainVec)
        {
            double[] rveCoordinatesSmallStrainVec = ((transformationMatrix *(new Vector( smallStrainVec))).Data);
            return rveCoordinatesSmallStrainVec;
        }


        private (double[] , double[,] ) StressesAndConstitutiveMatrixTransformation( double[] trueStressVec, double[,] constitutiveMat)
        {            
            var transformedTrueStressVec = (transformationMatrix.Transpose() * (new Vector(trueStressVec))).Data;
            // TODO: CHECK matrix for corrupted data after multiplication
            var transformedConstitutiveMat = (transformationMatrix.Transpose() * (new Matrix2D(constitutiveMat)) * transformationMatrix).Data;

            return (transformedTrueStressVec, transformedConstitutiveMat);
        }

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    if (Math.Abs(Cijrs_prev[i, j] - constitutiveMatrix[i, j]) > 1e-10)
                        return true;

            return false;
        }



        #region IFiniteElementMaterial3D methodoi mia mia 

        public IMatrix2D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[3]); // TODOGerasimos arxiko constitutive mporei na upologizetai pio efkola
                return new Matrix2D( constitutiveMatrix.Data); // TODO: apla kratame to constitutive matrix san array[,] (alla matrix mporei na xrhimopoithei gia tis peristrofes)
            }
        }

        public double[] Stresses // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return  trueStressVec; }
        }

        public void SaveState()
        {
            //microAnalyzer.SaveMaterialState() //TODO this can be used as well but nessesary analyzers should be defined

            var linearSystems = new Dictionary<int, ILinearSystem>(); //TODO: these lines are not used 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            // TODO alternatively:
            //var linearSystems = new Dictionary<int, ILinearSystem>();
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //{
            //    linearSystems.Add(subdomain.ID, new SkylineLinearSystem(subdomain.ID, subdomain.Forces));// prosoxh sto Id twn subdomain
            //}
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdaterWithInitialConditions(model.Subdomains[0]) };          
            foreach (var subdomainUpdater in subdomainUpdaters)
            {
                subdomainUpdater.UpdateState();
            }
            // Alternatively
            //ILinearSystem[] linearSystemsArray = new[] { linearSystems[1] }; //TODO those depend on the number of model.SubdomainsDictionary.Values as well            
            //foreach (ILinearSystem subdomain in linearSystemsArray)
            //{
            //    subdomainUpdaters[linearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == subdomain.ID).Index].UpdateState();
            //}
        }

        public bool Modified
        {
            get { return modified; }
        }

        public void ResetModified()
        {
            modified = false;
        }

        public int ID
        {
            get { return 1000; }
        }


        #endregion
        // methodoi ews edw xrhsimopoiountai
        public void ClearState() 
        {
            // pithanws TODO 
        }
        public void ClearStresses()
        {
            // pithanws TODO 
        }
        public double[] Coordinates { get; set; }

        public double YoungModulus => throw new NotImplementedException(); // TODO: remove these from the interface.

        public double PoissonRatio => throw new NotImplementedException();


        #region transformation methods
        //TODO: implement and use methods for shell transformation
        private double[] transformTrueStressVec(double[] trueStressVec, double[] tangent1, double[] tangent2, double[] normal)
        {
            throw new NotImplementedException();
        }

        private double[,] TransformConstitutiveMatrix(double[,] constitutiveMat, double[] tangent1, double[] tangent2, double[] normal)
        {
            throw new NotImplementedException();
        }
        
        #endregion

        //todo delete when unnesessary
        #region Print methods for debug 
        private void PrintMethodsForDebug(double[][] KfpDq, double[][] f2_vectors, double[][] f3_vectors, double[][] KppDqVectors, double[][] f4_vectors, double[,] DqCondDq, double[,] d2W_dfdf, double[,] Cijrs)
        {
            string string0 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_kanonikh_gewmetria_2\REF2_10__000_renu_new_multiple_algorithms_check_develop_1GrSh_correct_coh_CHECK_integration\d2\";

            string string1 = String.Concat(string0, @"KfpDq_{0}.txt");

            for (int i2 = 0; i2 < KfpDq.GetLength(0); i2++)
            {
                string path = string.Format(string1, (i2 + 1).ToString());
                Vector data = new Vector(KfpDq[i2]);
                data.WriteToFile(path);
            }

            string string2 = String.Concat(string0, @"KffInvKfpDq_{0}.txt");



            for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
            {
                string path = string.Format(string2, (i2 + 1).ToString());
                Vector data = new Vector(f2_vectors[i2]);
                data.WriteToFile(path);
            }

            string string3 = String.Concat(string0, @"f3_vectors_{0}.txt");
            string string4 = String.Concat(string0, @"KppDqVectors_{0}.txt");
            string string5 = String.Concat(string0, @"f4_vectors_{0}.txt");

            for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
            {
                string path = string.Format(string3, (i2 + 1).ToString());
                Vector data = new Vector(f3_vectors[i2]);
                data.WriteToFile(path);

                path = string.Format(string4, (i2 + 1).ToString());
                data = new Vector(KppDqVectors[i2]);
                data.WriteToFile(path);

                path = string.Format(string5, (i2 + 1).ToString());
                data = new Vector(f4_vectors[i2]);
                data.WriteToFile(path);

            }

            PrintUtilities.WriteToFile(DqCondDq, String.Concat(string0, @"DqCondDq.txt"));
            PrintUtilities.WriteToFile(d2W_dfdf,  String.Concat(string0, @"d2W_dfdf.txt"));
            PrintUtilities.WriteToFile(Cijrs, String.Concat(string0, @"Cijrs.txt"));
        }
        #endregion


    }



}
