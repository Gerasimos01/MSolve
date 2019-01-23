using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.IGA.Tests;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.MultiscaleAnalysisMerge;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.SamplesConsole.SupportiveClasses;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ISAAR.MSolve.SamplesConsole
{
    class SeparatecodeCheckingClass3
    {
        public static void CheckHexaFirst()
        {
            var modelBuilder = new GrapheneReinforcedRVEBuilderExample3GrSh1RVEstifDegenAndLinearPeripheralHostEmbSearchDdm();
            var ModelAndNodes = modelBuilder.GetModelAndBoundaryNodes();
        }

        public static void CheckLargeModelSeparation()
        {
            var modelBuilder = new GrapheneReinforcedRVEBuilderExample35fe2boundstiffHostTestPostDataDdm(1);
            var ModelAndNodes = modelBuilder.GetModelAndBoundaryNodes();
        }
    }
}
