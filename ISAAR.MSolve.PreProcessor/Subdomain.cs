﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor.Interfaces;

namespace ISAAR.MSolve.PreProcessor
{
    public class Subdomain
    {
        //private readonly IList<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, Node> nodesDictionary = new Dictionary<int, Node>();
        private readonly Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly Dictionary<int, Dictionary<DOFType, int>> globalNodalDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private double[] forces;

        #region Properties
        //public IList<EmbeddedNode> EmbeddedNodes
        //{
        //    get { return embeddedNodes; }
        //}

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, Node> NodesDictionary
        {
            get { return nodesDictionary; }
        }

        public IList<Node> Nodes
        {
            get { return nodesDictionary.Values.ToList<Node>(); }
        }

        public Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary
        {
            get { return nodalDOFsDictionary; }
        }

        public Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary
        {
            get { return globalNodalDOFsDictionary; }
        }

        public double[] Forces
        {
            get { return forces; }
        }

        public bool MaterialsModified { get; set; }
        //public bool MaterialsModified
        //{
        //    get
        //    {
        //        bool modified = false;
        //        foreach (Element element in elementsDictionary.Values)
        //            if (element.ElementType.MaterialModified)
        //            {
        //                modified = true;
        //                break;
        //            }
        //        return modified;
        //    }
        //}

        public int ID { get; set; }
        public int TotalDOFs { get; set; }
        #endregion

        #region Data inteconnection routines
        public void EnumerateDOFs()
        {
            TotalDOFs = 0;
            Dictionary<int, List<DOFType>> nodalDOFTypesDictionary = new Dictionary<int, List<DOFType>>();
            foreach (Element element in elementsDictionary.Values)
            {
                for (int i = 0; i < element.Nodes.Count; i++)
                {
                    if (!nodalDOFTypesDictionary.ContainsKey(element.Nodes[i].ID))
                        nodalDOFTypesDictionary.Add(element.Nodes[i].ID, new List<DOFType>());
                    nodalDOFTypesDictionary[element.Nodes[i].ID].AddRange(element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
                }
            }

            foreach (Node node in nodesDictionary.Values)
            {
                //List<DOFType> dofTypes = new List<DOFType>();
                //foreach (Element element in node.ElementsDictionary.Values)
                //{
                //    if (elementsDictionary.ContainsKey(element.ID))
                //    {
                //        foreach (DOFType dof in element.ElementType.DOFTypes)
                //            dofTypes.Add(dof);
                //    }
                //}

                Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
                //foreach (DOFType dofType in dofTypes.Distinct<DOFType>())
                foreach (DOFType dofType in nodalDOFTypesDictionary[node.ID].Distinct<DOFType>())
                {
                    int dofID = 0;
                    foreach (DOFType constraint in node.Constraints)
                    {
                        if (constraint == dofType)
                        {
                            dofID = -1;
                            break;
                        }
                    }
                    //var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
                    ////if (node.EmbeddedInElement != null && node.EmbeddedInElement.ElementType.GetDOFTypes(null)
                    ////    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    ////    dofID = -1;
                    //if (embeddedNode != null && embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null)
                    //    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    //    dofID = -1;

                    if (dofID == 0)
                    {
                        dofID = TotalDOFs;
                        TotalDOFs++;
                    }
                    dofsDictionary.Add(dofType, dofID);
                }
                
                nodalDOFsDictionary.Add(node.ID, dofsDictionary); 
            }
            forces = new double[TotalDOFs];
        }

        public void AssignGlobalNodalDOFsFromModel(Dictionary<int, Dictionary<DOFType, int>> glodalDOFsDictionary)
        {
            foreach (int nodeID in nodalDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = nodalDOFsDictionary[nodeID];
                Dictionary<DOFType, int> globalDOFTypes = new Dictionary<DOFType,int>(dofTypes.Count);
                foreach (DOFType dofType in dofTypes.Keys) 
                    globalDOFTypes.Add(dofType, glodalDOFsDictionary[nodeID][dofType]);
                globalNodalDOFsDictionary.Add(nodeID, globalDOFTypes);
            }
        }

        public void BuildNodesDictionary()
        {
            List<int> nodeIDs = new List<int>();
            Dictionary<int, Node> nodes = new Dictionary<int, Node>();
            foreach (Element element in elementsDictionary.Values)
                foreach (Node node in element.Nodes)
                {
                    nodeIDs.Add(node.ID);
                    if (!nodes.ContainsKey(node.ID))
                        nodes.Add(node.ID, node);
                }

            nodeIDs = new List<int>(nodeIDs.Distinct<int>());
            nodeIDs.Sort();
            foreach (int nodeID in nodeIDs)
                nodesDictionary.Add(nodeID, nodes[nodeID]);

            //foreach (var e in modelEmbeddedNodes.Where(x => nodeIDs.IndexOf(x.Node.ID) >= 0))
            //    embeddedNodes.Add(e);
        }

        #endregion

        public int[] GetCornerNodes()
        {
            int nodex1y1z1 = -1;
            int nodex2y1z1 = -1;
            int nodex1y2z1 = -1;
            int nodex2y2z1 = -1;
            int nodex1y1z2 = -1;
            int nodex2y1z2 = -1;
            int nodex1y2z2 = -1;
            int nodex2y2z2 = -1;
            double x1 = Double.MaxValue;
            double x2 = Double.MinValue;
            double y1 = Double.MaxValue;
            double y2 = Double.MinValue;
            double z1 = Double.MaxValue;
            double z2 = Double.MinValue;

            foreach (var kv in nodesDictionary)
            {
                if (x1 > kv.Value.X) x1 = kv.Value.X;
                if (x2 < kv.Value.X) x2 = kv.Value.X;
                if (y1 > kv.Value.Y) y1 = kv.Value.Y;
                if (y2 < kv.Value.Y) y2 = kv.Value.Y;
                if (z1 > kv.Value.Z) z1 = kv.Value.Z;
                if (z2 < kv.Value.Z) z2 = kv.Value.Z;

                if (x1 == kv.Value.X && y1 == kv.Value.Y && z1 == kv.Value.Z) nodex1y1z1 = kv.Key;
                if (x2 == kv.Value.X && y1 == kv.Value.Y && z1 == kv.Value.Z) nodex2y1z1 = kv.Key;
                if (x1 == kv.Value.X && y2 == kv.Value.Y && z1 == kv.Value.Z) nodex1y2z1 = kv.Key;
                if (x2 == kv.Value.X && y2 == kv.Value.Y && z1 == kv.Value.Z) nodex2y2z1 = kv.Key;
                if (x1 == kv.Value.X && y1 == kv.Value.Y && z2 == kv.Value.Z) nodex1y1z2 = kv.Key;
                if (x2 == kv.Value.X && y1 == kv.Value.Y && z2 == kv.Value.Z) nodex2y1z2 = kv.Key;
                if (x1 == kv.Value.X && y2 == kv.Value.Y && z2 == kv.Value.Z) nodex1y2z2 = kv.Key;
                if (x2 == kv.Value.X && y2 == kv.Value.Y && z2 == kv.Value.Z) nodex2y2z2 = kv.Key;
            }

            return new[] { nodex1y1z1, nodex2y1z1, nodex1y2z1, nodex2y2z1, nodex1y1z2, nodex2y1z2, nodex1y2z2, nodex2y2z2 };
        }

        public double[] GetLocalVectorFromGlobal(Element element, double[] globalVector)
        {
            int localDOFs = 0;
            foreach (IList<DOFType> dofs in element.ElementType.DOFEnumerator.GetDOFTypes(element)) localDOFs += dofs.Count;
            double[] localVector = new double[localDOFs];

            int pos = 0;
            for (int i = 0; i < element.ElementType.DOFEnumerator.GetDOFTypes(element).Count; i++)
            {
                //Arxikh diatupwsi MSOLVE
                //Node node = element.Nodes[i];
                //Diorthosi diatupwsis MSOLVE
                Node node = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element)[i];
                //ews edw
                foreach (DOFType dofType in element.ElementType.DOFEnumerator.GetDOFTypes(element)[i])
                {
                    int dof = NodalDOFsDictionary[node.ID][dofType];
                    if (dof != -1) localVector[pos] = globalVector[dof];
                    pos++;
                }
            }
            return localVector;
        }

        public void AddLocalVectorToGlobal(Element element, double[] localVector, double[] globalVector)
        {
            int pos = 0;
            for (int i = 0; i < element.ElementType.DOFEnumerator.GetDOFTypes(element).Count; i++)
            {
                //Arxikh diatupwsi MSOLVE
                //Node node = element.Nodes[i];
                //Diorthosi diatupwsis MSOLVE
                Node node = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element)[i];
                //ews edw
                foreach (DOFType dofType in element.ElementType.DOFEnumerator.GetDOFTypes(element)[i])
                {
                    int dof = NodalDOFsDictionary[node.ID][dofType];
                    if (dof != -1) globalVector[dof] += localVector[pos];
                    pos++;
                }
            }
        }

        // prosthiki print
        int ekteleseis_counter = 0;
        string string1 = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output\U_sunol_{0}.txt";

        public IVector<double> GetRHSFromSolution(IVector<double> solution, IVector<double> dSolution)
        {
            ekteleseis_counter += 1;
            ////if (ekteleseis_counter == 1)
            ////{ solution.WriteToFile(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\U_sunol.txt"); }
            ////if (ekteleseis_counter == 2)
            ////{ solution.WriteToFile(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\U_sunol_2.txt"); }
            string counter_data = ekteleseis_counter.ToString();
            string path = string.Format(string1, counter_data);
            solution.WriteToFile(path);
            Vector<double> forces = new Vector<double>(TotalDOFs);
            int current_element = 0;
            foreach (Element element in elementsDictionary.Values)
            {
                double[] localSolution = GetLocalVectorFromGlobal(element, ((Vector<double>)solution).Data);
                double[] localdSolution = GetLocalVectorFromGlobal(element, ((Vector<double>)dSolution).Data);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified) 
                    element.Subdomain.MaterialsModified = true;
                double[] f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                AddLocalVectorToGlobal(element, f, forces.Data);
                current_element += 1;
                //if(current_element==12)
                //{
                //    forces.WriteToFile(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\Forces_sunol_mono_shell_hexa.txt");
                //}
            }
            ////if (ekteleseis_counter == 1)
            ////{ forces.WriteToFile(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\Forces_sunol.txt"); }
            ////if (ekteleseis_counter == 2)
            ////{ forces.WriteToFile(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\Forces_sunol_2.txt"); }
            ////if (ekteleseis_counter == 3)
            ////{ forces.WriteToFile(@"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_embeded_shell_gia_check_tou_rve_embedding_sto_MSolve\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1\Forces_sunol_converged_Msolve.txt"); }
            return forces;
        }

        public void SaveMaterialState()
        {
            foreach (Element element in elementsDictionary.Values) element.ElementType.SaveMaterialState();
        }

        public void ResetMaterialsModifiedProperty()
        {
            this.MaterialsModified = false;
            foreach (Element element in elementsDictionary.Values) element.ElementType.ResetMaterialModified();
        }

        public void ClearMaterialStresses()
        {
            foreach (Element element in elementsDictionary.Values) element.ElementType.ClearMaterialStresses();
        }
    }
}
