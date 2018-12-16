using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    public interface IdegenerateRVEbuilder:IRVEbuilder
    {
        Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model model);
    }
}
