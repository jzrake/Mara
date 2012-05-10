local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'

local Quiet = false
local RunArgs = {
   N           = 128,
   id          = "test",
   CFL         = 0.8,
   tmax        = 0.4,
   noplot      = false,
   eosfile     = "none", -- tabeos.h5
   fluid       = "euler",
   boundary    = "outflow",
   advance     = "rk4",
   riemann     = "hllc",
   godunov     = "weno-split",
   reconstruct = "weno5",
   eos         = "gamma-law",
   adgam       = 1.4
}

for k,v in pairs(cmdline.opts) do
   if type(RunArgs[k]) == 'number' then
      RunArgs[k] = tonumber(v)
   else
      RunArgs[k] = v
   end
end


local function HandleErrors(Status, attempt)
   return 0
end

local Euler1dProblems = {

   Shocktube1 = {
      Pl = { 1.000, 1.000, 0.000, 0.0, 0.0 },
      Pr = { 0.125, 0.100, 0.000, 0.0, 0.0 } },

   Shocktube2 = {
      Pl = { 1.000, 0.400,-2.000, 0.0, 0.0 },
      Pr = { 1.000, 0.400, 2.000, 0.0, 0.0 } },

   Shocktube3 = {
      Pl = { 1.0, 1e+3, 0.0, 0.0, 0.0 },
      Pr = { 1.0, 1e-2, 0.0, 0.0, 0.0 } },

   Shocktube4 = {
      Pl = { 1.0, 1e-2, 0.0, 0.0, 0.0 },
      Pr = { 1.0, 1e+2, 0.0, 0.0, 0.0 } },

   Shocktube5 = {
      Pl = { 5.99924, 460.894, 19.59750, 0.0, 0.0 },
      Pr = { 5.99924,  46.095, -6.19633, 0.0, 0.0 } },

   ContactWave = {
      Pl = { 1.0, 1.0, 0.0, 0.7, 0.2 },
      Pr = { 0.1, 1.0, 0.0, 0.7, 0.2 }
   },

   IsentropicPulse = {
      pinit = function(x,y,z)
                 local n = 4
                 local L = 1.0
                 local K = 0.1
                 local Gamma = 1.4
                 local rho_ref = 1.0
                 local pre_ref = K * rho_ref ^ Gamma
                 local cs_ref = (Gamma * pre_ref / rho_ref)^0.5

                 local function f(x)
                    return math.sin(n*math.pi*x/L)^2
                 end

                 local rho = rho_ref * (1.0 + f(x))
                 local pre = K * rho ^ Gamma
                 local cs = (Gamma * pre/rho)^0.5
                 local vx = 2 / (Gamma - 1) * (cs - cs_ref)
                 return { rho, pre, vx, 0, 0 }
              end,
      entropy = 0.1 -- set equal to K above
   }
}

local Rmhd1dProblems = {

   Shocktube1 = {
      Pl = { 1.000, 1.000, 0.000, 0.0, 0.0, 0.5, 1.0, 0.0 },
      Pr = { 0.125, 0.100, 0.000, 0.0, 0.0, 0.5,-1.0, 0.0 } },

   Shocktube2 = {
      Pl = { 1.080, 0.950, 0.400, 0.3, 0.2, 2.0, 0.3, 0.3 },
      Pr = { 1.000, 1.000,-0.450,-0.2, 0.2, 2.5,-0.7, 0.5 } },

   Shocktube3 = {
      Pl = { 1.000, 0.100, 0.999, 0.0, 0.0, 10.0, 0.7, 0.7 },
      Pr = { 1.000, 0.100,-0.999, 0.0, 0.0, 10.0,-0.7,-0.7 } },

   Shocktube4 = {
      Pl = { 1.000, 5.000, 0.000, 0.3, 0.4, 1.0, 6.0, 2.0 },
      Pr = { 0.900, 5.300, 0.000, 0.0, 0.0, 1.0, 5.0, 2.0 } },

   ContactWave = {
      Pl = { 10.0, 1.0, 0.0, 0.7, 0.2, 5.0, 1.0, 0.5 },
      Pr = {  1.0, 1.0, 0.0, 0.7, 0.2, 5.0, 1.0, 0.5 }
   },

   RotationalWave = {
      Pl = { 1, 1, 0.400000, -0.300000, 0.500000, 2.4, 1.00,-1.600000 },
      Pr = { 1, 1, 0.377347, -0.482389, 0.424190, 2.4,-0.10,-2.178213 }
   },
}

local function setup()
   local N = RunArgs.N
   local NumberOfGhosts = { rmhd=8, srhd=5, euler=5 }
   set_domain({0.0}, {1.0}, {N}, NumberOfGhosts[RunArgs.fluid], 3)
   set_fluid(RunArgs.fluid)
   set_boundary(RunArgs.boundary)
   set_advance(RunArgs.advance)
   set_riemann(RunArgs.riemann)
   set_godunov(RunArgs.godunov)
   set_reconstruct(RunArgs.reconstruct)
   set_eos(RunArgs.eos,RunArgs.adgam)
end

local function CompareWenoEuler()

   local problem = tests.MakeShocktubeProblem(Euler1dProblems.Shocktube1, {reverse=false})
   util.run_simulation(problem:get_pinit(), setup , RunArgs)

   local P = get_prim()

   if RunArgs.noplot ~= '1' then
      util.plot{rho=P.rho, pre=P.pre, vy=P.vy, vz=P.vz}
   end
end

local function IsentopicConvergenceRate()

   RunArgs.boundary = "periodic"

   util.run_simulation(Euler1dProblems.IsentropicPulse.pinit, setup , RunArgs)

   local P = get_prim()
   if RunArgs.noplot ~= '1' then
      util.plot{rho=P.rho}
   end

   if false then
      local fout = io.open("rho.dat", "w")
      for i=0,RunArgs.N-1 do
	 fout:write(P.rho[i].."\n")
      end
   end
end

local function CompareEosRmhd()

   RunArgs.fluid   = "rmhd"
   RunArgs.riemann = "hlld"
   RunArgs.advance = "single"
   RunArgs.godunov = "plm-muscl"

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

--CompareEosRmhd()
CompareWenoEuler()
--IsentopicConvergenceRate()
