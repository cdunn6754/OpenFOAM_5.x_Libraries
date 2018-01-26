10-23-17:

I am really implementing the devolatilization stuff now (only primary for the moment) and I am making the assumption that the rate we calculate for volatile mass transfer is based on the daf particle weight and not the normal particle weight. 

Really we solve for the rate of change of the massfraction and then multiply it over the time step. So once we have the total change in mass fraction it needs to be multiplied by a mass so that we can return a "mass to transfer to carrier phase" from the calculate function. 

Just wanted to document an assumption there. I feel like I have read it somewhere in the documentation for PCCL and it makes sense since every other mass fraction or mass based quantity is based on daf mass not particle mass itself. 

