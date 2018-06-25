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


namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public class Microstructure2 //: IFiniteElementMaterial3D
    {
        private Model model { get; set; }
        //private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private Dictionary<int, Node> boundaryNodes { get; set; }
        private IRVEbuilder rveBuilder;
        private NewtonRaphsonNonLinearAnalyzer microAnalyzer;

        public Microstructure2(IRVEbuilder rveBuilder)
        {
            this.rveBuilder = rveBuilder;
            Tuple<Model, Dictionary<int, Node>> modelAndBoundaryNodes = this.rveBuilder.GetModelAndBoundaryNodes();
            this.model = modelAndBoundaryNodes.Item1;
            this.boundaryNodes = modelAndBoundaryNodes.Item2;
            this.model.ConnectDataStructures();
        }

        public object Clone()
        {
            return new Microstructure(rveBuilder);
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

            //if (matrices_not_initialized) TODOGerasimos
            //{ this.InitializeMatrices(); }            
            //D_tan_prev[k, j] = D_tan[k, j];

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            var solver = new SolverSkyline(linearSystems[1]);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 1; // microAnalyzer onomazetai o childAnalyzer
            microAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            // epivolh metakinhsewn ston analyzer pou molis dhmiourghsame
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                double[,] Dq_nodal = new double[9, 3];
                Dq_nodal[0, +0] = model.NodesDictionary[boundaryNode.ID].X; // h kai katedtheian boundaryNode.X 
                Dq_nodal[1, +1] = model.NodesDictionary[boundaryNode.ID].Y;
                Dq_nodal[2, +2] = model.NodesDictionary[boundaryNode.ID].Z;
                Dq_nodal[3, +0] = model.NodesDictionary[boundaryNode.ID].Y;
                Dq_nodal[4, +1] = model.NodesDictionary[boundaryNode.ID].Z;
                Dq_nodal[5, +2] = model.NodesDictionary[boundaryNode.ID].X;
                Dq_nodal[6, +0] = model.NodesDictionary[boundaryNode.ID].Z;
                Dq_nodal[7, +1] = model.NodesDictionary[boundaryNode.ID].X;
                Dq_nodal[8, +2] = model.NodesDictionary[boundaryNode.ID].Y;

                double[] u_prescr_xyz = new double[3];

                for (int i1 = 0; i1 < 3; i1++)
                {
                    for (int j1 = 0; j1 < 9; j1++)
                    {
                        u_prescr_xyz[i1] += Dq_nodal[j1, i1] * DefGradVec[j1]; //einai sunolikh //TODO
                    }
                }


            }
            // epivolh metakinhsewn ews edw
            microAnalyzer.SetMaxIterations = 100;
            microAnalyzer.SetIterationsForMatrixRebuild = 1;

            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, microAnalyzer, linearSystems);
            //childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 }); // MS na xthsimopoihsoume LOGFactories gia th dunamh isws?

            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            //TODOGerasimos to Solve o parent analyzer to pernaei ston NRNLAnalyzer o opoios kanei thn epanalhptikh diadikasia kai emeis tha tou exoume
            // kanei comment out to SaveMaterialState and update solution kai tha to kanoume Manually sto save state. (mono to save state oxi kai to u Add du)

            //TODOGerasimos
            double[] FPK_vec = new double[9];
            //TODO: upologismos FPK_vec apo oloklhrwma 
            double[,] DefGradMat = new double[3, 3] { { DefGradVec[0], DefGradVec[3], DefGradVec[6] }, { DefGradVec[7], DefGradVec[1], DefGradVec[4] }, { DefGradVec[5], DefGradVec[8], DefGradVec[2] } };
            double[,] FPK_mat = new double[3, 3] { { FPK_vec[0], FPK_vec[3], FPK_vec[6] }, { FPK_vec[7], FPK_vec[1], FPK_vec[4] }, { FPK_vec[5], FPK_vec[8], FPK_vec[2] } };
            double[,] SPK_mat = transformFPKtoSPK(DefGradMat, FPK_mat);

            //na elegxthei h parapanw anadiataxh kai o pollaplasiasmos
            double[,] d2W_dfdf = new double[9, 9];
            //TODO: upologismos d2W_dfdf apo oloklhrwma 



            //this.modified = CheckIfConstitutiveMatrixChanged();


        }

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



    }



}
