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



namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class Microstructure3Develop : IFiniteElementMaterial3D
    {
        private Model model { get; set; }
        //private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private Dictionary<int, Node> boundaryNodes { get; set; }
        private IRVEbuilder rveBuilder;
        //private NewtonRaphsonNonLinearAnalyzer microAnalyzer;
        private double volume;
        Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain;
        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements;
        private IScaleTransitions scaleTransitions = new DefGradVec3DScaleTransition();


        // aparaithta gia to implementation tou IFiniteElementMaterial3D
        IMatrix2D constitutiveMatrix;
        private double[] SPK_vec;
        private bool modified; // opws sto MohrCoulomb gia to modified

        private double[,] Cijrs_prev;
        private bool matrices_not_initialized = true;
        private double tol;
        public void InitializeMatrices()
        {
            Cijrs_prev = new double[6, 6];
            matrices_not_initialized = false;
            tol = Math.Pow(10, -19);
        }


        //double[] Stresses { get; }
        //IMatrix2D ConstitutiveMatrix { get; } TODOGerasimos

        public Microstructure3Develop(IRVEbuilder rveBuilder)
        {
            this.rveBuilder = rveBuilder;
            Tuple<Model, Dictionary<int, Node>,double> modelAndBoundaryNodes = this.rveBuilder.GetModelAndBoundaryNodes();
            this.model = modelAndBoundaryNodes.Item1;
            this.boundaryNodes = modelAndBoundaryNodes.Item2;
            this.volume = modelAndBoundaryNodes.Item3;
            this.model.ConnectDataStructures();
            this.InitializeFreeAndPrescribedDofsInitialDisplacementVectors();
        }

        private void InitializeFreeAndPrescribedDofsInitialDisplacementVectors()
        {
            uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, Vector>();
            foreach(Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                uInitialFreeDOFDisplacementsPerSubdomain.Add(subdomain.ID, new Vector(subdomain.TotalDOFs));// prosoxh sto Id twn subdomain
            }
            double[,] DGtr = new double[3, 3] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] DefGradVec = new double[9] { DGtr[0, 0], DGtr[1, 1], DGtr[2, 2], DGtr[1, 0], DGtr[2, 1], DGtr[0, 2], DGtr[2, 0], DGtr[0, 1], DGtr[1, 2], };
            initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                DefGradVec, initialConvergedBoundaryDisplacements);
            }

            throw new NotImplementedException();
        }

        public object Clone()
        {
            return new Microstructure2(rveBuilder);
        }

        public Dictionary<int, Node> BoundaryNodesDictionary
        {
            get { return boundaryNodes; }
        }
        public IList<Node> BoundaryNodes
        {
            get { return boundaryNodes.Values.ToList<Node>(); }
        }

        public void UpdateMaterial(double[] DefGradVec)
        {

            if (matrices_not_initialized) 
            { this.InitializeMatrices(); } // mporei na mpei kai ston constructor to initialization

            for (int i1 = 0; i1 < 6; i1++)
            {
                for (int j1 = 0; j1 < 6; j1++)
                {Cijrs_prev[i1, j1] = constitutiveMatrix[i1, j1];}
            }

            #region Creation of nessesary analyzers for NRNLAnalyzer

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            // TODO alternatively:
            //var linearSystems = new Dictionary<int, ILinearSystem>();
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //{
            //    linearSystems.Add(subdomain.ID, new SkylineLinearSystem(subdomain.ID, subdomain.Forces));// prosoxh sto Id twn subdomain
            //}

            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            var solver = new SolverSkyline(linearSystems[1]);
            var linearSystemsArray = new[] { linearSystems[1] }; //TODO those depend on the number of model.SubdomainsDictionary.Values as well
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdaterWithInitialConditions(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
            Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));
            #endregion

            #region Rve prescribed Dofs total DIsplacement Dictionary Creation (nessesary for NRNLAnalyzer)
            // epivolh metakinhsewn ston analyzer pou molis dhmiourghsame --> tha ginetai sto dictionary
            Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                scaleTransitions.ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(boundaryNode,
                DefGradVec, totalPrescribedBoundaryDisplacements);
            }
            #endregion


            #region Creation of Microstructure analyzer (NRNLdevelop temporarilly). 
            var increments = 1; // microAnalyzer onomazetai o childAnalyzer
            NewtonRaphsonNonLinearAnalyzerDevelopCopy microAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelopCopy(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
                boundaryNodes, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, equivalentContributionsAssemblers);
            microAnalyzer.SetMaxIterations = 100;
            microAnalyzer.SetIterationsForMatrixRebuild = 1;
            #endregion

            #region solution and update of free and prescribed converged displacements vectors;
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, microAnalyzer, linearSystems);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            uInitialFreeDOFDisplacementsPerSubdomain = microAnalyzer.GetConvergedSolutionVectorOfFreeDofs();
            initialConvergedBoundaryDisplacements = totalPrescribedBoundaryDisplacements;
            #endregion


            #region INTEGRATION stresses 
            double[] FppReactionVector = SubdomainCalculations.CalculateFppReactionsVector(Subdomain subdomain, Subdomain2 subdomain2, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, IVector solution, IVector dSolution)

            double[] FPK_vec = new double[9];
            //TODO: upologismos FPK_vec apo oloklhrwma 
            double[,] DefGradMat = new double[3, 3] { { DefGradVec[0], DefGradVec[3], DefGradVec[6] }, { DefGradVec[7], DefGradVec[1], DefGradVec[4] }, { DefGradVec[5], DefGradVec[8], DefGradVec[2] } };
            double[,] FPK_mat = new double[3, 3] { { FPK_vec[0], FPK_vec[3], FPK_vec[6] }, { FPK_vec[7], FPK_vec[1], FPK_vec[4] }, { FPK_vec[5], FPK_vec[8], FPK_vec[2] } };
            double[,] SPK_mat = transformFPKtoSPK(DefGradMat, FPK_mat);
            SPK_vec = new double[6] { SPK_mat[0,0], SPK_mat[1,1], SPK_mat[2,2], SPK_mat[0,1], SPK_mat[1,2], SPK_mat[0,2] };
            //na elegxthei h parapanw anadiataxh kai o pollaplasiasmos
            #endregion

            #region INTEGRATION constitutive Matrix
            //TODOGERASIMOS edw thelei NEWtonRaphson.BuildMatrices kai sigoura to solver. Initialize kapou edw
            microAnalyzer.BuildMatrices();
            double[][] KfpDq = SubdomainCalculations.CalculateKfreeprescribedDqMultiplications(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);

            // to BUildmatirces to exoume krathsei panw exw apo th sunarthsh tou f2
            double[][] f2_vectors = SubdomainCalculations.CalculateKffinverseKfpDq(KfpDq, model, elementProvider, scaleTransitions, boundaryNodes, solver, linearSystems);            
            double[][] f3_vectors = SubdomainCalculations.CalculateKpfKffinverseKfpDq(f2_vectors, model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);
            double[][] KppDqVectors = SubdomainCalculations.CalculateKppDqMultiplications(model.Subdomains[0], elementProvider, scaleTransitions, boundaryNodes);
            double[][] f4_vectors = SubdomainCalculations.SubtractConsecutiveVectors(KppDqVectors, f3_vectors);
            double[,] DqCondDq = SubdomainCalculations.CalculateDqCondDq(f4_vectors, scaleTransitions, boundaryNodes);

            double[,] d2W_dfdf = new double[DqCondDq.GetLength(0), DqCondDq.GetLength(1)];
            for (int i1 = 0; i1 < DqCondDq.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < DqCondDq.GetLength(1); i2++)
                {
                    d2W_dfdf[i1, i2] = (1 / volume) * DqCondDq[i1, i2];
                }
            }
            #endregion


            double[,] Cinpk = Transform_d2Wdfdf_to_Cijrs(d2W_dfdf, SPK_mat, DefGradMat); // to onomazoume Cinpk epeidh einai to 9x9 kai to diakrinoume etsi apo to Cijrs 6x6
            // transformation se 6x6 se 2 vhmata
            double[,] Cijrs_columns = new double[9, 6];
            for (int i1 = 0; i1 < 9; i1++)
            {
                Cijrs_columns[i1, 0] = Cinpk[i1, 0];
                Cijrs_columns[i1, 1] = Cinpk[i1, 1];
                Cijrs_columns[i1, 2] = Cinpk[i1, 2];
                Cijrs_columns[i1, 3] = 0.5 * (Cinpk[i1, 3] + Cinpk[i1, 7]);
                Cijrs_columns[i1, 4] = 0.5 * (Cinpk[i1, 4] + Cinpk[i1, 8]);
                Cijrs_columns[i1, 5] = 0.5 * (Cinpk[i1, 5] + Cinpk[i1, 6]);
            }

            double[,] Cijrs = new double[6, 6];

            for (int j1 = 0; j1 < 6; j1++)
            {
                Cijrs[0, j1] = Cijrs_columns[0, j1];
                Cijrs[1, j1] = Cijrs_columns[1, j1];
                Cijrs[2, j1] = Cijrs_columns[2, j1];
                Cijrs[3, j1] = 0.5 * (Cijrs_columns[3, j1] + Cijrs_columns[7, j1]);
                Cijrs[4, j1] = 0.5 * (Cijrs_columns[4, j1] + Cijrs_columns[8, j1]);
                Cijrs[5, j1] = 0.5 * (Cijrs_columns[5, j1] + Cijrs_columns[6, j1]);
            }

            constitutiveMatrix= new Matrix2D(Cijrs);

            this.modified = CheckIfConstitutiveMatrixChanged(); 
        }


        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (Math.Abs(Cijrs_prev[i, j] - constitutiveMatrix[i, j]) > 1e-10)
                        return true;

            return false;
        }



        #region IFiniteElementMaterial3D methodoi mia mia 

        public IMatrix2D ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[9] { 1,1,1,0,0,0,0,0,0 }); // TODOGerasimos arxiko constitutive mporei na upologizetai pio efkola
                return constitutiveMatrix;
            }
        }

        public double[] Stresses // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return SPK_vec; }
        }

        public void SaveState()
        {
            //microAnalyzer.SaveMaterialState() //TODO this can be used as well but nessesary analyzers should be defined

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
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







        #region transformation methods
        private double[,] transformFPKtoSPK(double[,] DefGradMat, double[,] FPK_mat)
        {
            IMatrix2D DefGradMat2D = new Matrix2D(DefGradMat);
            int linearsystemID = 1;
            SkylineLinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[3] { FPK_mat[0, 0], FPK_mat[1, 0], FPK_mat[2, 0] });
            //ILinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[3] { FPK_mat[0, 0], FPK_mat[1, 0], FPK_mat[2, 0] });
            var solver = new SolverSkyline(linearSystem);

            // BuildMatrices();
            linearSystem.Matrix = DefGradMat2D;

            //solver.Initialize();
            solver.Initialize(); // dld factorize

            // me thn parakatw commented out entolh anathetoume to linearSystem.RHS se ena vector sto opoio mporoumen na anaferthoume na kanoume copyTo klp.
            //Vector linearSystemRHS = ((Vector)linearSystem.RHS); // opws fainetai sth diadikasia pou kanei o NRNLAnalyzer sthn clculateInternalRHS sto telos tou loop

            double[,] SPK_Mat = new double[3, 3];
            Vector solution = new Vector(new double[3]);

            for (int j1 = 0; j1 < 3; j1++)
            {
                Vector RHS = new Vector(new double[3] { FPK_mat[0, j1], FPK_mat[1, j1], FPK_mat[2, j1] });
                SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
                k.Solve(RHS, solution);
                for (int i1 = 0; i1 < 3; i1++)
                {
                    SPK_Mat[i1, j1] = solution[i1];
                }
            }

            return SPK_Mat;
        }

        private double[,] Transform_d2Wdfdf_to_Cijrs(double[,] Aijkl, double[,] SPK, double[,] F)
        {
            int[,] i_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };
            int[,] k_seira = { { 1, 2, 3 }, { 3, 1, 2 }, { 2, 3, 1 } };

            double[,] Cinpk = new double[9, 9];

            double[,] F__F__ = new double[9, 9];
            //[F(:, 1) * F(1,:), F(:, 2) * F(1,:), F(:, 3) * F(1,:);
            //F(:, 1) * F(2,:),F(:, 2) * F(2,:),F(:, 3) * F(2,:);
            //F(:, 1) * F(3,:),F(:, 2) * F(3,:),F(:, 3) * F(3,:)];

            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int j1 = 0; j1 < 3; j1++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            F__F__[3 * i1 + k, 3 * j1 + l] = F[k, j1] * F[i1, l];
                        }
                    }

                }
            }

            //TODO upologismos F_inv (F__F__^(-1))

            double[,] F_inv = new double[9, 9];
            double[,] eye9 = new double[9, 9];
            for (int i1 = 0; i1 < 9; i1++)
            { eye9[i1, i1] = 1; }
            Vector solution = new Vector(new double[9]);

            SkylineMatrix2D F__F__Mat = new SkylineMatrix2D(F__F__);
            int linearsystemID = 1;
            SkylineLinearSystem linearSystem = new SkylineLinearSystem(linearsystemID, new double[9]);
            var solver = new SolverSkyline(linearSystem);
            // BuildMatrices();
            linearSystem.Matrix = F__F__Mat;
            //solver.Initialize();
            solver.Initialize(); // dld factorize

            for (int j1 = 0; j1 < 9; j1++)
            {
                Vector RHS = new Vector(new double[9] { eye9[0, j1], eye9[1, j1], eye9[2, j1], eye9[3, j1], eye9[4, j1], eye9[5, j1], eye9[6, j1], eye9[7, j1], eye9[8, j1] });
                SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
                k.Solve(RHS, solution);
                for (int i1 = 0; i1 < 9; i1++)
                {
                    F_inv[i1, j1] = solution[i1];
                }
            }


            for (int i1 = 0; i1 < 3; i1++)
            {
                for (int k1 = 0; k1 < 3; k1++)
                {
                    double[] A_j_l = new double[9] { Aijkl[3 * (i_seira[i1-1, 0] - 1) + i1, 3 * (k_seira[k1-1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 1] - 1) + i1, 3 * (k_seira[k1-1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 2] - 1) + i1, 3 * (k_seira[k1-1, 0] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 0] - 1) + i1, 3 * (k_seira[k1-1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 1] - 1) + i1, 3 * (k_seira[k1-1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 2] - 1) + i1, 3 * (k_seira[k1-1, 1] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 0] - 1) + i1, 3 * (k_seira[k1-1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 1] - 1) + i1, 3 * (k_seira[k1-1, 2] - 1) + k1],
                    Aijkl[3 * (i_seira[i1-1, 2] - 1) + i1, 3 * (k_seira[k1-1, 2] - 1) + k1]};

                    double[] sec_term = new double[9] { -SPK[i1 - 1, k1 - 1], 0, 0, 0, -SPK[i1 - 1, k1 - 1], 0, 0, 0, -SPK[i1 - 1, k1 - 1] };

                    Matrix2D F_invMat = new Matrix2D(F_inv);
                    Vector A_j_lVec = new Vector(A_j_l);
                    Vector sec_termVec = new Vector(sec_term);
                    Vector C_np_ = F_invMat * (new Vector(A_j_lVec + sec_termVec));

                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[0];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[1];
                    Cinpk[3 * (i_seira[i1, 0] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[2];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[3];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[4];
                    Cinpk[3 * (i_seira[i1, 1] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[5];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[0, k1] - 1) + 1] = C_np_[6];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[1, k1] - 1) + 2] = C_np_[7];
                    Cinpk[3 * (i_seira[i1, 2] - 1) + i1, 3 * (k_seira[2, k1] - 1) + 3] = C_np_[8];

                }
            }

            return Cinpk;
        }
        #endregion 


    }



}
