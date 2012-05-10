local json = require 'json'
local host = require 'host'
local util = require 'util'
local tests = require 'tests'

local Quiet = false
local RunArgs = {
   N           = 16,
   id          = "test",
   CFL         = 0.4,
   tmax        = 6.0,
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
   interactive = true
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





function ExplosionRmhd(x,y,z)
   local r2 = x*x + y*y
   if r2 < 0.025 then
      return { 1.000, 1.0, 0, 0, 0, 4, 0, 0 }
   else
      return { 0.125, 0.1, 0, 0, 0, 4, 0, 0 }
   end
end


local function setup()
   local N = RunArgs.N
   local NumberOfConserved = { rmhd=8, srhd=5, euler=5 }
--   set_domain({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, {N,N,N}, NumberOfGhosts[RunArgs.fluid], 7
  set_domain({-0.5, -0.5}, {0.5, 0.5}, {N,N}, NumberOfConserved[RunArgs.fluid], 3)
   set_fluid(RunArgs.fluid)
   set_eos(RunArgs.eos,RunArgs.adgam)
   set_boundary(RunArgs.boundary)
   set_advance(RunArgs.advance)
   set_riemann(RunArgs.riemann)
   set_godunov(RunArgs.godunov)
end

util.run_simulation(tests.KelvinHelmoltz:get_pinit(), setup , RunArgs)
