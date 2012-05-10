local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'

local RunArgs = {
   N           = 64,
   dim         = 1,
   id          = "test",
   ic          = "Shocktube1", -- name of test problem
   CFL         = 0.5,
   tmax        = 0.2,
   noplot      = false,
   eosfile     = "none", -- tabeos.h5
   fluid       = "euler",
   boundary    = "outflow",
   advance     = "rk4",
   riemann     = "hllc",
   godunov     = "weno-split",
   reconstruct = "weno5",
   eos         = "gamma-law",
   adgam       = 1.4,
   vis         = "default",
   quiet       = false,
   problem     = "defaultshocktube",
   interactive = false
}

for k,v in pairs(cmdline.opts) do
   if type(RunArgs[k]) == 'number' then
      RunArgs[k] = tonumber(v)
   else
      RunArgs[k] = v
   end
end

if RunArgs.dim ~= 1 and RunArgs.vis == "default" then
   RunArgs.interactive = true
elseif RunArgs.dim == 1 and RunArgs.vis == "default" then
   RunArgs.interactive = false
elseif RunArgs.vis == "false" then
   RunArgs.interactive = false
elseif RunArgs.vis == "true" then
   RunArgs.interactive = true
else print("vis must equal 'false' or 'true'")
end


local function HandleErrors(Status, attempt)
   return 0
end

local function setup()
   local N = RunArgs.N
   local NumberOfConserved = { rmhd=8, srhd=5, euler=5 }

   if RunArgs.dim == 1 then
      set_domain({0.0}, {1.0}, {N},
                 NumberOfConserved[RunArgs.fluid], 3)
   elseif RunArgs.dim == 2 then
      set_domain({-0.5, -0.5}, {0.5, 0.5}, {N,N},
                 NumberOfConserved[RunArgs.fluid], 3)
   elseif RunArgs.dim == 3 then
      set_domain({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, {N,N,N},
                 NumberOfConserved[RunArgs.fluid], 3)
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




----------------------------------------------------------
local function DensityWaveConvergenceRate()
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
----------------------------------------------------------
local function IsentopicConvergenceRate()
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
----------------------------------------------------------
local function CompareEosRmhd()
   RunArgs.fluid   = "rmhd"
   RunArgs.riemann = "hlld"
   RunArgs.advance = "single"
   RunArgs.godunov = "plm-muscl"
   RunArgs.boundary = "outflow"

   if RunArgs.eosfile ~= "none" then
      local tabeos = require 'tabeos'
      tabeos.MakeNeutronStarUnits()
      tabeos.LoadMicroPh(RunArgs.eosfile)
   end

   local Dl = 2.0
   local Dr = 1.0
   local Tl = 8.0 -- MeV
   local Tr = 2.0 -- MeV
   local prel = eos.Pressure(Dl, Tl)
   local prer = eos.Pressure(Dr, Tr)

   local ShocktubeEos = {
      Pl = { Dl, prel, 0.000, 0.0, 0.0, 0.0, 0.0, 0.0 },
      Pr = { Dr, prer, 0.000, 0.0, 0.0, 0.0, 0.0, 0.0 } }

   local function pinit(x,y,z)
      if x < 0.5 then
         return ShocktubeEos.Pl
      else
         return ShocktubeEos.Pr
      end
   end

   util.run_simulation(pinit, setup , RunArgs)

   local P = get_prim()
   if RunArgs.noplot ~= '1' then
      util.plot{rho=P.rho, pre=P.pre*1000}
   end
end

----------------------------------------------------------
local function RmhdExplosion()
   RunArgs.fluid   = "rmhd"
   RunArgs.riemann = "hlld"
   RunArgs.advance = "single"
   RunArgs.godunov = "plm-muscl"
   RunArgs.boundary = "outflow"
   util.run_simulation(tests.Explosion:get_pinit(), setup , RunArgs)

   local P = get_prim()
   if RunArgs.dim == 1 then
      util.plot{rho=P.rho, pre=P.pre, vx=P.vx, vy=P.vy, vz=P.vz}
   end
end

----------------------------------------------------------
local function VanillaShocktube()

   local problem = tests.MakeShocktubeProblem(tests[RunArgs.ic], {reverse=false})
   util.run_simulation(problem:get_pinit(), setup , RunArgs)

   local P = get_prim()
   if RunArgs.dim == 1 then
      util.plot{rho=P.rho, pre=P.pre, vx=P.vx, vy=P.vy, vz=P.vz}
   end
end
----------------------------------------------------------
local function VanillaExplosion()
   util.run_simulation(tests.Explosion:get_pinit(), setup , RunArgs)

   local P = get_prim()
   if RunArgs.dim == 1 then
      util.plot{rho=P.rho, pre=P.pre, vx=P.vx, vy=P.vy, vz=P.vz}
   end
end
----------------------------------------------------------
local function VanillaKelvinHelmoltz()
   util.run_simulation(tests.KelvinHelmoltz:get_pinit(), setup , RunArgs)

   local P = get_prim()
   if RunArgs.dim == 1 then
      util.plot{rho=P.rho, pre=P.pre, vx=P.vx, vy=P.vy, vz=P.vz}
   end
end
----------------------------------------------------------
local function VanillaDensityWave()
   util.run_simulation(tests.DensityWave:get_pinit(), setup , RunArgs)

   local P = get_prim()
   if RunArgs.dim == 1 then
      util.plot{rho=P.rho, pre=P.pre, vx=P.vx, vy=P.vy, vz=P.vz}
   end
end
----------------------------------------------------------
local function VanillaIsentropicPulse()
   util.run_simulation(tests.IsentropicPulse:get_pinit(), setup , RunArgs)
   local P = get_prim()

   if RunArgs.dim == 1 then
      util.plot{rho=P.rho, pre=P.pre, vx=P.vx, vy=P.vy, vz=P.vz}
   end
end
----------------------------------------------------------

if RunArgs.problem == "defaultshocktube" then
   VanillaShocktube()
elseif RunArgs.problem == "explosion" then
   VanillaExplosion()
elseif RunArgs.problem == "kelvinhelmoltz" then
   VanillaKelvinHelmoltz()
elseif RunArgs.problem == "densitywave" then
   VanillaDensityWave()
elseif RunArgs.problem == "isentropicpulse" then
   VanillaIsentropicPulse()
elseif RunArgs.problem == "IsentropicConvergence" then
   IsentopicConvergenceRate()
elseif RunArgs.problem == "DensityWaveConvergence" then
   DensityWaveConvergenceRate()
elseif RunArgs.problem == "CompareEosRmhd" then
   CompareEosRmhd()
elseif RunArgs.problem == "rmhdexplosion" then
   RmhdExplosion()
else
   print("Error: No such problem")
end
