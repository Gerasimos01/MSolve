﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Transfer;


//TODO: Transfer level sets and recalculate nodal enrichments in each process. Perhaps also identify which nodes are enriched 
//      with what. 
namespace ISAAR.MSolve.XFEM.Entities
{
    public class XModelMpi : ModelMpiBase<XModel>
    {
        //TODO: This does not guarantee that the model also uses the same elementFactory for the elements of this process's 
        //      subdomain.
        private readonly IXFiniteElementFactory elementFactory;

        public XModelMpi(ProcessDistribution processDistribution, Func<XModel> createModel, 
            IXFiniteElementFactory elementFactory) : base(processDistribution)
        {
            this.elementFactory = elementFactory;
            if (processDistribution.IsMasterProcess) this.model = createModel();
            else this.model = new XModel();
        }

        public IDomain2DBoundary Boundary => this.model.Boundary;

        public XSubdomain GetXSubdomain(int subdomainID)
        {
            procs.CheckProcessMatchesSubdomainUnlessMaster(subdomainID);
            return model.Subdomains[subdomainID];
        }

        public override void ScatterSubdomains()
        {
            // Serialize the data of each subdomain
            XSubdomainDto[] serializedSubdomains = null;
            if (procs.IsMasterProcess)
            {
                serializedSubdomains = new XSubdomainDto[procs.Communicator.Size];

                for (int p = 0; p < procs.Communicator.Size; ++p)
                {
                    if (p == procs.MasterProcess) serializedSubdomains[p] = XSubdomainDto.CreateEmpty();
                    else
                    {
                        XSubdomain subdomain = model.Subdomains[procs.GetSubdomainIdOfProcess(p)];
                        serializedSubdomains[p] = XSubdomainDto.Serialize(subdomain, DofSerializer);
                    }
                }
            }

            // Scatter the serialized subdomain data from master process
            XSubdomainDto serializedSubdomain = procs.Communicator.Scatter(serializedSubdomains, procs.MasterProcess);

            // Deserialize and store the subdomain data in each process
            if (!procs.IsMasterProcess)
            {
                XSubdomain subdomain = serializedSubdomain.Deserialize(DofSerializer, elementFactory);
                model.Subdomains[subdomain.ID] = subdomain;
                subdomain.ConnectDataStructures();
            }
        }
    }
}
