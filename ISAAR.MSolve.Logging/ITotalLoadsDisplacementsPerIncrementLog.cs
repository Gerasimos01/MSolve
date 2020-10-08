using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging
{
    public interface ITotalLoadsDisplacementsPerIncrementLog
    {
        void LogTotalDataForIncrement(int incrementNumber, int currentIterationNumber, double errorNorm,
            IVectorView totalDisplacements, IVectorView totalInternalForces);

        void Initialize();

    }
}
