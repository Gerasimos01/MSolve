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
    class SeparateCodeCheckingClass4
    {
        public static void Check05bStressIntegrationObjeIntegration()
        {
            //PROELEFSI: SeparateCodeCheckingClass.Check05bStressIntegration
            //allages: tha xrhsimopoithei h nea microstructure me obje kapoia subdomainCalculations

            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            ElasticMaterial3D material1 = new ElasticMaterial3D()
            { YoungModulus = E_disp, PoissonRatio = ni_disp, };
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = SeparateCodeCheckingClass.Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = SeparateCodeCheckingClass.Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateMaterial(new StressStrainVectorContinuum3D(GLVec));
            material1.SaveState();
            double[] stressesCheck2 = material1.Stresses.Data;

            VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheck27Hexa();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            IContinuumMaterial3DDefGrad microstructure3 = new Microstructure3DevelopMultipleSubdomainsUseBaseSimuRandObj(homogeneousRveBuilder1,false,1);
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateMaterial(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses.Data;
            microstructure3.SaveState();
            microstructure3.UpdateMaterial(new double[9] { 1.20, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses.Data;
        }
    }
}
