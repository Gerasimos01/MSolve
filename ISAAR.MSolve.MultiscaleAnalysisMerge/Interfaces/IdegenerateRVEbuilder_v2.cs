using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    public interface IdegenerateRVEbuilder_v2 : IRVEbuilder_v2
    {
        Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model_v2 model);
    }
}
