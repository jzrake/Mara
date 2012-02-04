
local json = require 'json'
local host = require 'host'
local util = require 'util'


local RunArgs = {
   N       = 128,
   id      = "test",
   CFL     = 0.8,
   tmax    = 0.4,
   noplot  = false
}


local function HandleErrors(Status, attempt)
   return 0
end

-- *****************************************************************************
-- Main driver, operates between checkpoints and then returns
-- .............................................................................
local function RunSimulation(Status, Howlong)

   local t0 = Status.CurrentTime
   local attempt = 0

   while Status.CurrentTime - t0 < Howlong do

      if HandleErrors(Status, attempt) ~= 0 then
	 return 1
      end
      attempt = attempt + 1

      local dt = Status.Timestep
      local kzps, errors = advance(dt)

      if errors == 0 then
         driving.Advance(dt)
         driving.Resample()

         print(string.format("%05d(%d): t=%5.4f dt=%5.4e %3.2fkz/s %3.2fus/(z*Nq)",
                             Status.Iteration, attempt-1, Status.CurrentTime, dt,
			     kzps, 1e6/8/(1e3*kzps)))
	 io.flush()

	 attempt = 0
         Status.Timestep = get_timestep(RunArgs.CFL)
         Status.CurrentTime = Status.CurrentTime + Status.Timestep
         Status.Iteration = Status.Iteration + 1
      end
   end
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

local function InitSimulation(problem, setup)

   for k,v in pairs(cmdline.opts) do
      if type(RunArgs[k]) == 'number' then
	 RunArgs[k] = tonumber(v)
      else
	 RunArgs[k] = v
      end
   end

   print("runtime arguments:")
   for k,v in pairs(RunArgs) do
      print(k,v)
   end

   local function pinit(x,y,z)
      if x < 0.5 then
	 return problem.Pl
      else
	 return problem.Pr
      end
   end

   setup()
   init_prim(pinit)

   local Status = { }

   Status.CurrentTime = 0.0
   Status.Iteration   = 0
   Status.Checkpoint  = 0
   Status.Timestep    = 0.0

   local datadir = string.format("data/%s", RunArgs.id)
   os.execute(string.format("mkdir -p %s", datadir))
   os.execute(host.Filesystem(datadir))

   print_mara()
   return Status
end


local function setup_plm()
   local N = RunArgs.N
   set_domain({0.0}, {1.0}, {N}, 5, 2)
   set_fluid("euler")
   set_boundary("outflow")
   set_advance("single")
   set_riemann("hllc")
   set_godunov("plm-muscl", 2.0, 0)
   set_eos("gamma-law", 1.4)
end

local function setup_weno()
   local N = RunArgs.N
   set_domain({0.0}, {1.0}, {N}, 5, 3)
   set_fluid("euler")
   set_boundary("outflow")
   set_advance("rk3")
   set_godunov("weno-split")
   set_eos("gamma-law", 1.4)
end

local function setup_rmhd()
   local N = RunArgs.N
   set_domain({0.0}, {1.0}, {N}, 8, 2)
   set_fluid("rmhd")
   set_boundary("outflow")
   set_riemann("hlld")
   set_advance("single")
   set_godunov("plm-muscl", 2.0, 0)
   set_eos("gamma-law", 1.4)
end


local Status = InitSimulation(Euler1dProblems.Shocktube1, setup_plm)
RunSimulation(Status, RunArgs.tmax)
local P_plm = get_prim()

local Status = InitSimulation(Euler1dProblems.Shocktube1, setup_weno)
RunSimulation(Status, RunArgs.tmax)
local P_weno = get_prim()

if noplot == '1' then
   util.plot{weno=P_weno.rho, plm=P_plm.rho}
   util.plot{weno=P_weno.vz, plm=P_plm.vz}
end