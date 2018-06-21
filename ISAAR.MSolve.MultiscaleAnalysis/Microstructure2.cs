using System;
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

        public Microstructure2(IRVEbuilder rveBuilder )
        {
            this.rveBuilder = rveBuilder;
            Tuple < Model, Dictionary<int, Node> > modelAndBoundaryNodes=this.rveBuilder.GetModelAndBoundaryNodes();
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
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) }

            var increments = 1; // microAnalyzer onomazetai o childAnalyzer
            microAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            // epivolh metakinhsewn ston analyzer pou molis dhmiourghsame
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                double[,] Dq_nodal = new double[9, 3];
                Dq_nodal[0,  + 0] = model.NodesDictionary[boundaryNode.ID].X; // h kai katedtheian boundaryNode.X 
                Dq_nodal[1,  + 1] = model.NodesDictionary[boundaryNode.ID].Y;
                Dq_nodal[2,  + 2] = model.NodesDictionary[boundaryNode.ID].Z;
                Dq_nodal[3,  + 0] = model.NodesDictionary[boundaryNode.ID].Y;
                Dq_nodal[4,  + 1] = model.NodesDictionary[boundaryNode.ID].Z;
                Dq_nodal[5,  + 2] = model.NodesDictionary[boundaryNode.ID].X;
                Dq_nodal[6,  + 0] = model.NodesDictionary[boundaryNode.ID].Z;
                Dq_nodal[7,  + 1] = model.NodesDictionary[boundaryNode.ID].X;
                Dq_nodal[8,  + 2] = model.NodesDictionary[boundaryNode.ID].Y;

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
            double[,] DefGradMat = new double [3,3]{ { DefGradVec[0], DefGradVec[3], DefGradVec[6] }, { DefGradVec[7], DefGradVec[1], DefGradVec[4] }, { DefGradVec[5], DefGradVec[8], DefGradVec[2] } };

            //na elegxthei h parapanw anadiataxh kai o pollaplasiasmos
            double[,] d2W_dfdf= new double[9, 9];
            //TODO: upologismos d2W_dfdf apo oloklhrwma 



            //this.modified = CheckIfConstitutiveMatrixChanged();


        }



    }
}
