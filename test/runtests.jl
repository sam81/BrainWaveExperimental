t1 = time_ns()
#using BrainwaveExperimental
#using Base.Test

include("test_examples.jl")
# write your own tests here
#@test 1 == 1
t2 = time_ns() 
println(float(t2-t1)*1e9)
