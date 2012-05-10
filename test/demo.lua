local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'

local RunArgs = {
   N         = 64,
   dim       = 1,
   id        = "test",
   ic        = "Shocktube1", -- name of test problem
   CFL       = 0.6,
   tmax      = 0.2,
   noplot    = false,
   eosfile   = "none", -- tabeos.h5
   fluid     = "euler",
   boundary  = "outflow",
   advance   = "rk4",
   riemann   = "hllc",
   godunov   = "weno-split",

   -- --------------------------------------------------------------------------
   -- flags sent to config_solver
   -- --------------------------------------------------------------------------
   fsplit    = "llf",  -- one of [llf, marq]       ... flux splitting mode
   extrap    = "weno5",-- one of [pcm, plm, weno5] ... reconstruction type
   theta     = 2.0,    -- must be [0,2]            ... theta value for PLM/minmod
   IS        = "js96", -- one of [js96, b08, sz10] ... smoothness indicator
   sz10A     = 50.0,   -- should be in [0,100]     ... used by sz10 (see weno.c)

   -- --------------------------------------------------------------------------
   -- flags for particular initial conditions setups
   -- --------------------------------------------------------------------------
   eos       = "gamma-law",
   adgam     = 1.4,
   quiet     = false,
   problem   = "shocktube",
   plotvars  = "rho,pre,vx,vy,vz",
   angle     = "{1,0,0}",
}
util.parse_args(RunArgs)
tests.RunArgs = RunArgs


local function HandleErrors(Status, attempt)
   return 0
end

local function plot_prim()
   if not RunArgs.noplot then
      local P = get_prim()
      local pltdict =  { }
      for k,v in pairs(util.string_split(RunArgs.plotvars, ",")) do
	 pltdict[v] = P[v]
      end
      if RunArgs.dim == 1 then
	 util.plot(pltdict)
      end
   end
end

local function cfg_mara()
   local N = RunArgs.N
   local NumberOfConserved = { rmhd=8, srhd=5, euler=5 }
   local Nq = NumberOfConserved[RunArgs.fluid]

   if RunArgs.dim == 1 then
      set_domain({0}, {1}, {N}, Nq, 3)
   elseif RunArgs.dim == 2 then
      set_domain({0,0}, {1,1}, {N,N}, Nq, 3)
   elseif RunArgs.dim == 3 then
      set_domain({0,0,0}, {1,1,1}, {N,N,N}, Nq, 3)
   else
      error("Invalid Dimension")
   end
   
   if type(RunArgs.boundary) == "string" then
      set_boundary(RunArgs.boundary)
   elseif type(RunArgs.boundary) == "table" then
      set_boundary(unpack(RunArgs.boundary))
   end
   config_solver(RunArgs)
   set_fluid(RunArgs.fluid)
   set_advance(RunArgs.advance)
   set_riemann(RunArgs.riemann)
   set_godunov(RunArgs.godunov)
   set_eos(RunArgs.eos,RunArgs.adgam)
end


local ProblemList = { }

function ProblemList.DensityWaveConvergenceRate()
   local outf = io.open("densitywave.dat", "w")
   RunArgs.boundary = "periodic"

   if RunArgs.dim == 1 then
      local res_values = { 64, 128, 256, 512, 1024 }
   else
      local res_values = { 16, 32, 64, 128}
   end

   for run_num,N in pairs(res_values) do
      RunArgs.N = N
      local problem = tests.DensityWave
      problem.eps = 3.2e-1

      problem.velocity = { 0.1, 0.0, 0.0 }

      local status = util.run_simulation(problem:get_pinit(), cfg_mara, RunArgs)
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


function ProblemList.IsentopicConvergenceRate()
   local outf = io.open("isentropic.dat", "w")
   RunArgs.boundary = "periodic"

   if RunArgs.dim == 1 then
      local res_values = { 64, 128, 256, 512, 1024 }
   else
      local res_values = { 16, 32, 64, 128}
   end

   for run_num,N in pairs(res_values) do
      RunArgs.N = N
      local problem = tests.IsentropicPulse
      problem.mode = 2
      util.run_simulation(problem:get_pinit(), cfg_mara, RunArgs)

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



function ProblemList.RmhdExplosion()
   RunArgs.fluid   = "rmhd"
   RunArgs.riemann = "hlld"
   RunArgs.advance = "single"
   RunArgs.godunov = "plm-muscl"
   util.run_simulation(tests.Explosion:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.ImplosionProblem()
   RunArgs.fluid = "euler"
   RunArgs.boundary = {"reflect2d", 2, 3}
   util.run_simulation(tests[RunArgs.ic]:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaShocktube()
   local problem = tests[RunArgs.ic]
   if not problem then
      error("ic="..RunArgs.ic.." does not exist")
   end
   util.run_simulation(tests[RunArgs.ic]:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaExplosion()
   util.run_simulation(tests.Explosion:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaKelvinHelmholtz()
   RunArgs.dim = 2
   util.run_simulation(tests.KelvinHelmholtz:get_pinit(), cfg_mara , RunArgs)
end
function ProblemList.VanillaDensityWave()
   util.run_simulation(tests.DensityWave:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end
function ProblemList.VanillaIsentropicPulse()
   util.run_simulation(tests.IsentropicPulse:get_pinit(), cfg_mara , RunArgs)
   if RunArgs.dim == 1 then plot_prim() end
end


-- -----------------------------------------------------------------------------
-- Some shortcuts
-- -----------------------------------------------------------------------------
ProblemList["mhdexpl"] = ProblemList["RmhdExplosion"]
ProblemList["implode"] = ProblemList["ImplosionProblem"]
ProblemList["denswave"] = ProblemList["VanillaDensityWave"]
ProblemList["isenwave"] = ProblemList["VanillaIsentropicPulse"]
ProblemList["kh"] = ProblemList["VanillaKelvinHelmholtz"]
ProblemList["shocktube"] = ProblemList["VanillaShocktube"]

local prob_to_run = ProblemList[RunArgs.problem]
if not prob_to_run then
   error("problem="..RunArgs.problem.." not found")
else
   prob_to_run()
end
