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

-- *****************************************************************************
-- Main driver, operates between checkpoints and then returns
-- .............................................................................
function RunSimulation(CFL, EndTime)

   visual.open_window()

   local dt = 0.0
   local Iteration = 0
   local CurrentTime = 0.0
   local piter = true
   local start = os.clock()

   print_mara()
   init_prim(Explosion)--KelvinHelmoltz)

--   local problem = tests.DensityWave
--   problem.velocity = { 1.0, 1.0, 0 }
--   problem.mode = { 1.0, 1.0, 0 }
--   init_prim(tests.DensityWave:get_pinit(0.0))

   while CurrentTime < EndTime do

      local prim = get_prim()
      local draw_array = prim.rho[':,8,:']
      visual.draw_texture(draw_array)

--      local ppmname = string.format("images/output-%04d.ppm", Iteration)
--      write_ppm(ppmname, prim["rho"])
--      os.execute("convert -rotate 90 " .. ppmname .. " " .. string.gsub(ppmname, ".ppm", ".png"))
--      os.execute("rm " .. ppmname)

      local kzps = advance(dt)

      if piter then
         print(string.format("%05d: t=%5.4f, dt=%5.4e %3.2fkz/s",
                             Iteration, CurrentTime, dt, kzps))
      end

      dt = get_timestep(CFL)
      CurrentTime = CurrentTime + dt
      Iteration = Iteration + 1

      collectgarbage()
   end

   print(string.format("run took %3.2f seconds", os.clock() - start))
   return get_prim()
end


function Explosion(x,y,z)
   local r2 = x*x + y*y + z*z
   if r2 < 0.005 then
      return { 1.000, 1.0, 0, 0, 0 }
   else
      return { 0.125, 0.1, 0, 0, 0 }
   end
end


function KelvinHelmoltz(x,y,z)

   if math.abs(y) > 0.25 then
      rho = 1.0
      vx = -0.5
   else
      rho = 2.0
      vx =  0.5
   end

   vx = 0.02*(math.random() - 0.5) + vx
   vy = 0.02*(math.random() - 0.5)

   return { rho, 2.5, vx, vy, 0.0 }
end


function ExplosionRmhd(x,y,z)
   local r2 = x*x + y*y
   if r2 < 0.025 then
      return { 1.000, 1.0, 0, 0, 0, 4, 0, 0 }
   else
      return { 0.125, 0.1, 0, 0, 0, 4, 0, 0 }
   end
end

local N = RunArgs.N
local NumberOfGhosts = { rmhd=8, srhd=5, euler=5 }

set_domain({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, {N,N,N}, NumberOfGhosts[RunArgs.fluid], 7)
set_fluid(RunArgs.fluid)
set_eos(RunArgs.eos,RunArgs.adgam)
set_boundary(RunArgs.boundary)
set_advance(RunArgs.advance)
set_riemann(RunArgs.riemann)
set_godunov(RunArgs.godunov)
set_advance("rk4")
set_godunov("weno-split")

RunSimulation(RunArgs.CFL, RunArgs.tmax)
