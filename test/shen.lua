

json = require 'json'
host = require 'host'

local shen_file = host.ShenFile
tab = load_shen(shen_file, 0.08, {1e12, 1e16}, {0.1, 200.0})

--f = io.open("shen_08.json", "r")
--contents = f:read("*all")
--tab = json.decode(contents)

set_eos("shen", tab)

test_shen()
