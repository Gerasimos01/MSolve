using System;
using System.Collections.Generic;
//using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.FEM
{
    public static class SubdomainCalculationsCopyDebugMethods
    {
        public static double[][] CalculateKfreeprescribedDqMultiplications(Subdomain subdomain, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            //ElementStructuralStiffnessProvider s = new ElementStructuralStiffnessProvider();
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

            

            double[][] KfpDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
            for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
            {
                KfpDqVectors[j1] = new double[subdomain.TotalDOFs]; // h allliws subdomain.Forces.GetLength(0)
            }

            // TODO: should encapsulate DOF logic into a separate entity that will manage things if embedded or not (should return element matrix and globaldofs correspondence list
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            //SkylineMatrix2D K = new SkylineMatrix2D(GlobalMatrixAssemblerSkyline.CalculateRowIndex(subdomain, nodalDOFsDictionary)); // commented out den kanoume assembly afto
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);

            foreach (Element element in subdomain.ElementsDictionary.Values) // TODOGerasimos edw mporei na xrhsimopoihthei to dictionary twn eleement pou exoun fp nodes
            {
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    foreach (DOFType dofTypeRow in elementDOFTypes[i])
                    {
                        int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (dofRow != -1) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                        {                    // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                //{
                                //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                //    if (dofColumn != -1)
                                //    {
                                //        int height = dofRow - dofColumn;
                                //        if (height >= 0)
                                //            K.Data[K.RowIndex[dofRow] + height] += ElementK[iElementMatrixRow, iElementMatrixColumn];
                                //    }
                                //    iElementMatrixColumn++;
                                //}
                                int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
                                if (boundaryNodes.ContainsKey(nodeColumn.ID)) 
                                {
                                    double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 
                                    for (int j1 = 0; j1 < nodalDofsNumber; j1++)                  //scaleTransitions.MicroToMacroTransition ti tha xrhsimopoihsei apo afta analoga pws einai implemented
                                    {
                                       element_Kfp_triplette[j1]= ElementK[iElementMatrixRow, iElementMatrixColumn+j1];
                                    }

                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kfp_triplette);
                                    for (int j2 = 0; j2 < contribution.GetLength(0); j2++)
                                    {
                                        KfpDqVectors[j2][dofRow] += contribution[j2]; // TODO diorthothike
                                    }

                                }
                                iElementMatrixColumn += nodalDofsNumber;

                            }
                        }
                        iElementMatrixRow++;
                    }
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return KfpDqVectors;
        }

        public static double[][] CalculateKffinverseKfpDq(double[][] KfpDq, Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, SolverSkyline solver, Dictionary<int, ILinearSystem> linearSystems)
        {
            double[][] f2_vectors = new double[KfpDq.GetLength(0)][];
            for (int i1 = 0; i1 < KfpDq.GetLength(0); i1++)
            {
                //double[] FirstRHS = KfpDq[0];
                //int FirstRHSdimension = FirstRHS.GetLength;


                SkylineLinearSystem Kff_linearSystem = new SkylineLinearSystem(1, new double[model.Subdomains[0].TotalDOFs]);
                var K_ffsolver = new SolverSkyline(Kff_linearSystem);
                // BuildMatrices(); exei proigithei prin thn CalculateKfreeprescribedDqMultiplications klp
                Kff_linearSystem.Matrix = linearSystems[1].Matrix; // pairnoume dld to matrix apo ekei pou to topothetei o StaticAnalyzer otan kanei InitializeMatrices();
                //solver.Initialize();
                solver.Initialize(); // prin to factorize periexei if opote den tha kanei thn idia douleia duo fores

                Vector Kff_solution = new Vector(new double[model.Subdomains[0].TotalDOFs]);
                Vector K_ffRHS = new Vector(KfpDq[i1]);
                SkylineMatrix2D k_temp = ((SkylineMatrix2D)Kff_linearSystem.Matrix); // opws sto solverskyline.cs sthn Solve()
                k_temp.Solve(K_ffRHS, Kff_solution);
                f2_vectors[i1] = Kff_solution.Data;
            }

            return f2_vectors;
        }

        public static double[][] CalculateKpfKffinverseKfpDq(double[][] f2_vectors, Subdomain subdomain, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;
            double[][] f3_vectors = new double[f2_vectors.GetLength(0)][];
            for (int i1 = 0; i1 < f2_vectors.GetLength(0); i1++)
            {
                f3_vectors[i1] = new double[scaleTransitions.PrescribedDofsPerNode() * boundaryNodes.Count];
            }
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                {
                                    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                    if (dofColumn != -1)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                                    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
                                        
                                        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
                                        {
                                            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow+i1, iElementMatrixColumn]*f2_vectors[i2][dofColumn];
                                        }

                                    }
                                    iElementMatrixColumn++;
                                }
                            }

                        }
                    }
                    iElementMatrixRow+= elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return f3_vectors;

        }

        public static double[][] CalculateKppDqMultiplications(Subdomain subdomain, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary; //TODOGerasimos mallon commented out afto den xreiazetai

            double[][] KppDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);
            for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
            {
                KppDqVectors[j1] = new double[boundaryNodesOrder.Count*scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)
            }


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                //
                                //foreach (DOFType dofTypeColumn in elementDOFTypes[j])
                                //{
                                //    int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                //    if (dofColumn != -1)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                                //    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok

                                //        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
                                //        {
                                //            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow + i1, iElementMatrixColumn] * f2_vectors[i2][dofColumn]; /////
                                //        }

                                //    }
                                //    iElementMatrixColumn++;
                                //}
                                //
                                int nodalDofsNumber= elementDOFTypes[j].Count;
                                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                                {
                                    double[] element_Kpp_triplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                                    for (int j2 = 0; j2 < scaleTransitions.PrescribedDofsPerNode(); j2++)
                                    {
                                        element_Kpp_triplette[j2] = ElementK[iElementMatrixRow+i1, iElementMatrixColumn + j2]; // mallon iElementMatrixRow + i1
                                    }
                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kpp_triplette);
                                    for (int j1 = 0; j1 < contribution.GetLength(0); j1++)
                                    {
                                        KppDqVectors[j1][dofrow_p] += contribution[j1];
                                    }
                                }
                                iElementMatrixColumn += nodalDofsNumber;



                            }

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return KppDqVectors;
        }

        public static double[][] SubtractConsecutiveVectors(double[][] KppDqVectors,double[][] f3_vectors)
        {
            double[][] f4_vectors = new double[KppDqVectors.GetLength(0)][];
            for (int j1 = 0; j1 < f4_vectors.GetLength(0); j1++)
            {
                f4_vectors[j1] = new double[KppDqVectors[0].GetLength(0)];
                for (int j2 = 0; j2 < f4_vectors[0].GetLength(0); j2++)
                {
                    f4_vectors[j1][j2] = KppDqVectors[j1][j2] - f3_vectors[j1][j2];
                }
            }

            return f4_vectors;

        }

        public static double[,] CalculateDqCondDq(double[][] f4_vectors, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            double[,] DqCondDq = new double[scaleTransitions.MacroscaleVariableDimension(), scaleTransitions.MacroscaleVariableDimension()];

            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);

            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                for (int i1 = 0; i1 < f4_vectors.GetLength(0); i1++)
                {
                    double[] f4DataTriplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                    for (int i2 = 0; i2 < scaleTransitions.PrescribedDofsPerNode(); i2++)
                    {
                        f4DataTriplette[i2] = f4_vectors[i1][scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[boundaryNode.ID] - 1) + i2];
                    }
                    double[] contribution = scaleTransitions.MicroToMacroTransition(boundaryNode, f4DataTriplette);
                    for (int i3 = 0; i3 < scaleTransitions.MacroscaleVariableDimension(); i3++)
                    {
                       DqCondDq[i3, i1] += contribution[i3];
                    }
                }

            }
            return DqCondDq;
        }

        public static Dictionary<int,int> GetNodesOrderInDictionary(Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, int> boundaryNodesOrder = new Dictionary<int, int>();
            int order = 1;
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                boundaryNodesOrder.Add(boundaryNode.ID, order);
                order += 1;
            }
            return boundaryNodesOrder;

        }

        public static double[] CalculateFppReactionsVector(Subdomain subdomain,  IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, IVector solution, IVector dSolution,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            //TODOGerasimos: 1) Subdomain2 einai h upo kataskevh subdomain.cs ths Marias gia na mporoume na anaferthoume sthn methodo ths CalculateElementNodalDisplacements(..,..). 
            // Otan parei telikh morfh tha taftizetai me thn Subdomain.cs
            // 2)IVector solution, IVector dSolution EINAI AFTA ME TA OPOIA kaloume thn GetRHSFromSolution sthn 213 tou NRNLAnalyzer
            double[] FppReactionVector ;
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);
            FppReactionVector = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement; 
                var elStart = DateTime.Now;
                //IMatrix2D ElementK = elementProvider.Matrix(element);
                var localSolution = subdomain.GetLocalVectorFromGlobal(element, solution);
                subdomain.ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                double[] localdSolution = subdomain.GetLocalVectorFromGlobal(element, dSolution);
                double[] f = element.ElementType.CalculateForces(element, localSolution, localdSolution);

                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            FppReactionVector[dofrow_p] += f[iElementMatrixRow + i1];
                         
                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return FppReactionVector;

        }

        //public static double[] CalculateFppReactionsVectorCopy(Subdomain subdomain, Subdomain2 subdomain2, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, IVector solution, IVector dSolution)
        //{
        //    //TODOGerasimos: 1) Subdomain2 einai h upo kataskevh subdomain.cs ths Marias gia na mporoume na anaferthoume sthn methodo ths CalculateElementNodalDisplacements(..,..). 
        //    // Otan parei telikh morfh tha taftizetai me thn Subdomain.cs
        //    // 2)IVector solution, IVector dSolution EINAI AFTA ME TA OPOIA kaloume thn GetRHSFromSolution sthn 213 tou NRNLAnalyzer
        //    double[] FppReactionVector;
        //    Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);
        //    FppReactionVector = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)


        //    var times = new Dictionary<string, TimeSpan>();
        //    var totalStart = DateTime.Now;
        //    times.Add("rowIndexCalculation", DateTime.Now - totalStart);
        //    times.Add("element", TimeSpan.Zero);
        //    times.Add("addition", TimeSpan.Zero);
        //    foreach (Element element in subdomain.ElementsDictionary.Values)
        //    {
        //        var isEmbeddedElement = element.ElementType is IEmbeddedElement;
        //        var elStart = DateTime.Now;
        //        //IMatrix2D ElementK = elementProvider.Matrix(element);
        //        double[] localSolution = subdomain2.CalculateElementNodalDisplacements(element, solution);
        //        double[] localdSolution = subdomain2.CalculateElementNodalDisplacements(element, dSolution);
        //        double[] f = element.ElementType.CalculateForces(element, localSolution, localdSolution);

        //        times["element"] += DateTime.Now - elStart;

        //        elStart = DateTime.Now;
        //        var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
        //        var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
        //        int iElementMatrixRow = 0;
        //        for (int i = 0; i < elementDOFTypes.Count; i++)
        //        {
        //            Node nodeRow = matrixAssemblyNodes[i];
        //            if (boundaryNodes.ContainsKey(nodeRow.ID))
        //            {
        //                for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
        //                {
        //                    int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
        //                    FppReactionVector[dofrow_p] += f[iElementMatrixRow + i1];

        //                }
        //            }
        //            iElementMatrixRow += elementDOFTypes[i].Count;
        //        }
        //        times["addition"] += DateTime.Now - elStart;
        //    }
        //    var totalTime = DateTime.Now - totalStart;

        //    return FppReactionVector;

        //}


        public static double[] CalculateDqFpp(double[] FppReactionVector, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            double[] DqFpp = new double[scaleTransitions.MacroscaleVariableDimension()];

            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);

            foreach (Node boundaryNode in boundaryNodes.Values)
            {

                double[] FppDataTriplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                for (int i2 = 0; i2 < scaleTransitions.PrescribedDofsPerNode(); i2++)
                {
                    FppDataTriplette[i2] = FppReactionVector[scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[boundaryNode.ID] - 1) + i2];
                }
                double[] contribution = scaleTransitions.MicroToMacroTransition(boundaryNode, FppDataTriplette);
                for (int i3 = 0; i3 < scaleTransitions.MacroscaleVariableDimension(); i3++)
                {
                    DqFpp[i3] += contribution[i3];
                }


            }

            return DqFpp;

        }
    }
}
