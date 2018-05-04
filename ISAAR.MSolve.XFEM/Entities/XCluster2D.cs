﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// Inter-subdomain communication 
    /// </summary>
    class XCluster2D
    {
        private readonly List<XSubdomain2D> subdomains;

        public XCluster2D()
        {
            this.subdomains = new List<XSubdomain2D>();
        }

        public XClusterDofOrderer DofOrderer { get; private set; }
        public IReadOnlyList<XSubdomain2D> Subdomains { get { return subdomains; } }


        public void AddSubdomain(XSubdomain2D subdomain)
        {
            subdomains.Add(subdomain);
        }

        public void OrderDofs(Model2D model)
        {
            DofOrderer = XClusterDofOrderer.CreateNodeMajor(model, this);
        }
    }
}
