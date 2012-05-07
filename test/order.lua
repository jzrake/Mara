
local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'


local RunArgs = {
   N       = 128,
   id      = "test",
   CFL     = 0.6,
   tmax    = 0.2,
   noplot  = false,
   quiet   = false
}
util.parse_args(RunArgs)




local function IsentopicConvergenceRate()
   local outf = io.open("isentropic.dat", "w")
   local res_values = { 64, 128, 256, 512, 1024 }

   for run_num,N in pairs(res_values) do
      local function setup()
	 set_domain({0.0}, {1.0}, {N}, 5, 7)
	 set_fluid("euler")
	 set_boundary("periodic")
	 set_riemann("hllc")
--	 set_advance("single")
--	 set_godunov("plm-muscl", 2.0, 0)
	 set_advance("rk4")
	 set_godunov("weno-riemann")
--	 set_godunov("weno-split")
	 set_eos("gamma-law", 1.4)
      end

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
      local function setup()
	 set_domain({0.0}, {1.0}, {N}, 5, 7)
	 set_fluid("euler")
	 set_boundary("periodic")
	 set_riemann("hllc")
--	 set_advance("single")
--	 set_godunov("plm-muscl", 2.0, 0)
	 set_advance("rk4")
--	 set_godunov("weno-riemann")
	 set_godunov("weno-split")
	 set_eos("gamma-law", 1.4)
      end

      local problem = tests.DensityWave
      problem.eps = 3.2e-1
      local status = util.run_simulation(problem:get_pinit(), setup, RunArgs)
      
      local P_comp = get_prim()

      init_prim(problem:get_pinit(status.CurrentTime))
      local P_true = get_prim()

--      util.plot{rho_true=P_true.rho, rho_comp=P_comp.rho}
--      os.exit()

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
      local function setup()
	 set_domain({0.0, 0.0}, {1.0, 1.0}, {N, N}, 5, 3)
	 set_fluid("euler")
	 set_boundary("periodic")
	 set_riemann("hllc")
--	 set_advance("single")
--	 set_godunov("plm-muscl", 2.0, 0)
	 set_advance("rk4")
--	 set_godunov("weno-riemann")
	 set_godunov("weno-split")
	 set_eos("gamma-law", 1.4)
      end

      local problem = tests.DensityWave
      problem.eps = 3.2e-1
      local status = util.run_simulation(problem:get_pinit(), setup, RunArgs)
      
      local P_comp = get_prim()

      init_prim(problem:get_pinit(status.CurrentTime))
      local P_true = get_prim()

--      util.plot{rho_true=P_true.rho, rho_comp=P_comp.rho}
--      os.exit()

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
