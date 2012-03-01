

json = require 'json'
host = require 'host'

local shen_file = host.ShenFile
tab = load_shen(shen_file, 0.08, {1e12, 1e14}, {10.0, 100.0})

for k,v in pairs(tab) do
--   print (k,v:shape('array'))
--   print(v)
end

set_eos("shen", tab)
test_shen()


