﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Creates a elastic matrix rve
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class HomogeneousRVEBuilderLinear : IRVEbuilder
    {        
        private Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp ;
        private rveMatrixParameters mp;
        private grapheneSheetParameters gp;
        private string renumbering_vector_path;

        public HomogeneousRVEBuilderLinear()
        {
            //TODOGerasimos
            // this.renumbering_vector_path=renumbering_vector_path,
            // this.subdiscr1=subdiscr1
        }
        public IRVEbuilder Clone(int a) => new HomogeneousRVEBuilderLinear();

        public Tuple<Model, Dictionary<int, Node>,double> GetModelAndBoundaryNodes()
        {
           return Reference2RVEExample10_000withRenumbering_mono_hexa();
        }

        public ISolver GetAppropriateSolver(Model model)
        {
            return (new SkylineSolver.Builder()).BuildSolver(model);
        }


        public Tuple<Model, Dictionary<int, Node>,double> Reference2RVEExample10_000withRenumbering_mono_hexa()
        {
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain( 1 ));

            Dictionary<int, Node> boundaryNodes= new Dictionary<int, Node>();
            // COPY APO: Reference2RVEExample100_000withRenumbering_mono_hexa
            double[,] Dq = new double[1, 1];
            //Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            //rveMatrixParameters mp;
            //grapheneSheetParameters gp;
            renumbering_vector_path = "..\\..\\..\\RveTemplates\\Input\\RveHomogeneous\\REF_new_total_numbering27.txt";
            
            string Fxk_p_komvoi_rve_path = @"C:\Users\turbo-x\Desktop\notes_elegxoi\REFERENCE_fe2_diafora_check\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF2_10_000_renu_new_multiple_algorithms_check_stress_27hexa\Fxk_p_komvoi_rve.txt";
            int subdiscr1 = 1;
            int discr1 = 3;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 3;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;

            mpgp = FEMMeshBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 0; // 0 gra sheets afou exoume mono hexa
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];


            FEMMeshBuilder.LinearHexaElementsOnlyRVEwithRenumbering_forMS(model, mp, Dq, renumbering_vector_path, boundaryNodes);
            double volume = mp.L01 * mp.L02 * mp.L03;

            

            return new Tuple<Model, Dictionary<int, Node>,double>(model, boundaryNodes,volume);
        }

        


    }
    //HomogeneousRVEBuilderCheck27HexaLinear
    //TODOGerasimos gia na ta krataei mesa kai na kanei build model oses fores tou zhththei
    // omoiws na ginei kai to RVE me graphene sheets 
    // string renumbering_vector_path; 
    // int subdiscr1;
}
