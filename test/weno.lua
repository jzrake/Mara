


local matrix = require 'matrix'
local sep = "---------------------------------"
local function clip(x)
   if math.abs(x) < 1e-12 then
      return 0.0
   else
      return x
   end
end


-- Flux Jacobian as given by http://www.microcfd.com/download/pdf/aiaa2001-2609.pdf
-- .............................................................................
local function Jacobian1(U, n) -- n must be a unit vector

   local gm = 1.4
   local g1 = gm - 1.0
   local g2 = gm - 2.0

   local D = U[0] -- inputs are Mara convention: {D,E,px,py,pz}
   local E = U[1]
   local u = U[2]/D
   local v = U[3]/D
   local w = U[4]/D

   local nx, ny, nz = n[1], n[2], n[3]
   local vn = u*nx + v*ny + w*nz
   local ek = 0.5*(u^2 + v^2 + w^2)
   local p0 = (E-D*ek)*g1 -- pressure
   local a2 = gm*p0/D     -- sound speed
   local h0 = a2 / g1 + ek

   A = { { 0, nx, ny, nz, 0 },
         { g1*ek*nx - u*vn,
           1*vn - g2*u*nx,
           u*ny - g1*v*nx,
           u*nz - g1*w*nx, g1*nx },

         { g1*ek*ny - v*vn,
           v*nx - g1*u*ny,
           1*vn - g2*v*ny,
           v*nz - g1*w*ny, g1*ny },

         { g1*ek*nz - w*vn,
           w*nx - g1*u*nz,
           w*ny - g1*v*nz,
           1*vn - g2*w*nz, g1*nz },

         { (g1*ek-h0)*vn,
           h0*nx - g1*u*vn,
           h0*ny - g1*v*vn,
           h0*nz - g1*w*vn, gm*vn } }

   return matrix(A) -- output is A := dF{D,px,py,pz,E}/d{D,px,py,pz,E}, Toro's convention
end

local function Jacobian2(U, n) -- n must be {1,0,0}
   local gm = 1.4
   local g1 = gm - 1.0
   local g2 = gm - 2.0

   local D = U[0] -- inputs are Mara convention: {D,E,px,py,pz}
   local E = U[1]
   local u = U[2]/D
   local v = U[3]/D
   local w = U[4]/D

   local ek = 0.5*(u^2 + v^2 + w^2)
   local p0 = (E-D*ek)*g1 -- pressure
   local a2 = gm*p0/D     -- sound speed
   local H = a2 / g1 + ek

   A = { { 0, 1, 0, 0, 0 },
         { g1*H - u^2 - a2, (3-gm)*u, -g1*v, -g1*w, g1 },
         { -u*v, v, u, 0, 0 },
         { -u*w, w, 0, u, 0 },
         { 0.5*u*(H*(gm-3) + g1*ek - a2), H - g1*u^2, -g1*u*v, -g1*u*w, gm*u } }

   return matrix(A) -- output is A := dF{D,px,py,pz,E}/d{D,px,py,pz,E}, Toro's convention
end


local function TestState(P, dim)
   print("\n\n")
   print(sep)
   print("Testing state on dim:", dim)
   print(sep)

   local n = { 0,0,0 }
   n[dim] = 1.0
   local U = fluid.PrimToCons(P)
   local Lv, Rv, lm = fluid.Eigensystem(P, dim)

   local A = Jacobian1(U, n)
   local Lv, Rv, lm = fluid.Eigensystem(P, dim)

--   print((Jacobian2(U, n) - A):replace(clip))

   local L = matrix(Lv)
   local R = matrix(Rv)

   print("\nR.L (expect I):")
   print(sep)
   print((R*L):replace(clip))

   print("\neigenvalues:")
   print(sep)
   print(lm)

   print("\nL.A.R (expect eigenvalues):")
   print(sep)
   print((L*A*R):replace(clip))
end



local P = { 1, 2, 1, 2, 3 }
TestState(P, 3)
