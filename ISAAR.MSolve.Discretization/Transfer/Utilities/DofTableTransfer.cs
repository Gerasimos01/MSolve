﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using MPI;

namespace ISAAR.MSolve.Discretization.Transfer.Utilities
{
    public class DofTableTransfer
    {
        private readonly IModel model;
        private readonly ProcessDistribution procs;

        private IEnumerable<ISubdomain> subdomainsToGatherFrom_master;
        private bool gatherFromThisSubdomain_slave;
        private DofTable dofsToSend_slave;

        public DofTableTransfer(IModel model, ProcessDistribution processDistribution)
        {
            this.model = model;
            this.procs = processDistribution;
        }

        public Dictionary<int, int> NumSubdomainDofs_master { get; }
        public Dictionary<int, DofTable> SubdomainDofOrderings_master { get; }

        public void DefineModelData_master(IEnumerable<ISubdomain> subdomainsToGatherFrom)
        {
            // Single out the subdomain of the master process, since no MPI transfers are required for it.
            int masterSubdomain = procs.GetSubdomainIdOfProcess(procs.MasterProcess);
            this.subdomainsToGatherFrom_master = subdomainsToGatherFrom.Where(sub => sub.ID != masterSubdomain);
        }

        public void DefineSubdomainData_slave(bool gatherFromThisSubdomain, DofTable dofsToGather)
        {
            this.gatherFromThisSubdomain_slave = gatherFromThisSubdomain;
            this.dofsToSend_slave = dofsToGather;
        }

        public void Transfer(int mpiTag)
        {
            var tableSerializer = new DofTableSerializer(model.DofSerializer);

            // Receive the corner dof ordering of each subdomain from the corresponding process
            if (procs.IsMasterProcess)
            {
                var serializedTables = new Dictionary<ISubdomain, int[]>();
                foreach (ISubdomain subdomain in subdomainsToGatherFrom_master)
                {
                    // Receive the corner dof ordering of each subdomain
                    int source = procs.GetProcessOfSubdomain(subdomain.ID);
                    //Console.WriteLine($"Process {comm.Rank} (master): Started receiving dof ordering from process {source}.");
                    serializedTables[subdomain] = MpiUtilities.ReceiveArray<int>(procs.Communicator, source, mpiTag);
                    //Console.WriteLine($"Process {comm.Rank} (master): Finished receiving dof ordering from process {source}.");
                }

                // After finishing with all comunications deserialize the received items. //TODO: This should be done concurrently with the transfers, by another thread.
                foreach (ISubdomain subdomain in subdomainsToGatherFrom_master)
                {
                    bool isModified = serializedTables.TryGetValue(subdomain, out int[] serializedTable);
                    if (isModified)
                    {
                        //Console.WriteLine($"Process {comm.Rank} (master): Started deserializing corner dof ordering of subdomain {subdomain.ID}.");
                        NumSubdomainDofs_master[subdomain.ID] = tableSerializer.CountEntriesOf(serializedTable);
                        SubdomainDofOrderings_master[subdomain.ID] = 
                            tableSerializer.Deserialize(serializedTable, model.GetNode);
                        //Console.WriteLine($"Process {comm.Rank} (master): Finished deserializing corner dof ordering of subdomain {subdomain.ID}.");
                        serializedTables.Remove(subdomain); // Free up some temporary memory.
                    }
                }
            }
            else
            {
                if (gatherFromThisSubdomain_slave)
                {
                    //Console.WriteLine($"Process {rank}: Started serializing dof ordering.");
                    int[] serializedTable = tableSerializer.Serialize(dofsToSend_slave);
                    //Console.WriteLine($"Process {rank}: Finished serializing dof ordering. Started sending it to master");
                    MpiUtilities.SendArray(procs.Communicator, serializedTable, procs.MasterProcess, mpiTag);
                    //Console.WriteLine($"Process {rank}: Finished sending corner dof ordering to master.");
                }
            }
        }
    }
}
