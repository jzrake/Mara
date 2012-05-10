local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'



local RunArgs = {
   N           = 64,
   dim         = 1,
   id          = "test",
   CFL         = 0.8,
   tmax        = 0.4,
   noplot      = false,
   eosfile     = "none", -- tabeos.h5
   fluid       = "euler",
   boundary    = "periodic",
   advance     = "rk4",
   riemann     = "hllc",
   godunov     = "weno-split",
   reconstruct = "weno5",
   eos         = "gamma-law",
   adgam       = 1.4,
   interactive = false,
   quiet       = false
}
util.parse_args(RunArgs)

local function setup()
   local N = RunArgs.N
   local NumberOfConserved = { rmhd=8, srhd=5, euler=5 }

   if RunArgs.dim == 1 then
      set_domain({0.0}, {1.0}, {N}, NumberOfConserved[RunArgs.fluid], 3)
   elseif RunArgs.dim == 2 then
      set_domain({-0.5, -0.5}, {0.5, 0.5}, {N,N}, NumberOfConserved[RunArgs.fluid], 3)
   elseif RunArgs.dim == 3 then
      set_domain({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, {N,N,N}, NumberOfConserved[RunArgs.fluid], 3)
   else print("Invalid Dimension")
   end

   set_fluid(RunArgs.fluid)
   set_boundary(RunArgs.boundary)
   set_advance(RunArgs.advance)
   set_riemann(RunArgs.riemann)
   set_godunov(RunArgs.godunov)
   set_reconstruct(RunArgs.reconstruct)
   set_eos(RunArgs.eos,RunArgs.adgam)
end


local function IsentopicConvergenceRate()
   local outf = io.open("isentropic.dat", "w")
   local res_values = { 64, 128, 256, 512, 1024 }

   for run_num,N in pairs(res_values) do
      local problem = tests.IsentropicPulse
      problem.mode = 2
      util.run_simulation(problem:get_pinit(), setup, RunArgs)
      
      local P = get_prim()
      local dS = problem:entropy(P.rho, P.pre) - problem.entropy_ref

      local L1 = 0.0
      for i=0,dS:size()-1 do
	 L1 = L1 + math.abs(dS[i]) / N
      end

      print("L1 = " .. math.log10(L1))
      outf:write(N .. " " .. L1 .. "\n")
   end
end

local function DensityWaveConvergenceRate()
   local outf = io.open("densitywave.dat", "w")
   local res_values = { 64, 128, 256, 512, 1024 }

   for run_num,N in pairs(res_values) do
      local problem = tests.DensityWave
      problem.eps = 3.2e-1
      local status = util.run_simulation(problem:get_pinit(), setup, RunArgs)
      
      local P_comp = get_prim()

      init_prim(problem:get_pinit(status.CurrentTime))
      local P_true = get_prim()

      local diff = P_comp.rho - P_true.rho

      local L1 = 0.0
      for i=0,diff:size()-1 do
	 L1 = L1 + math.abs(diff[i]) / N
      end

      print("L1 = " .. math.log10(L1))
      outf:write(N .. " " .. L1 .. "\n")
   end
end


local function DensityWaveConvergenceRate2d()
   local outf = io.open("densitywave.dat", "w")
   local res_values = { 16, 32, 64, 128 }

   for run_num,N in pairs(res_values) do
      local problem = tests.DensityWave
      problem.eps = 3.2e-1
      problem.velocity = { 0.1, 0.0, 0.0 }
      local status = util.run_simulation(problem:get_pinit(), setup, RunArgs)
      
      local P_comp = get_prim()

      init_prim(problem:get_pinit(status.CurrentTime))
      local P_true = get_prim()

      local diff = P_comp.rho - P_true.rho

      local L1 = 0.0
      for i=0,diff:size()-1 do
	 L1 = L1 + math.abs(diff[i]) / (N*N)
      end

      print("L1 = " .. math.log10(L1))
      outf:write(N .. " " .. L1 .. "\n")
   end
end


--IsentopicConvergenceRate()
--DensityWaveConvergenceRate()
DensityWaveConvergenceRate2d()
