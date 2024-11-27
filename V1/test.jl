using Plots
using UnROOT

# Plot random data with labels
scatter(rand(10), rand(10), title="Random Data", xlabel="X", ylabel="Y")

println("Hello, World!")

# Directory of the ROOT file
dir = "rootfiles/run_01.root"

# Open the ROOT file
file = ROOTFile(dir)

# List keys in the ROOT file
println("Keys in ROOT file: ", keys(file))

# Access the "Delphes" tree
tree = file["Delphes"]

# List properties of the tree
println("Properties of 'Delphes' tree: ", propertynames(tree["Jet/Jet.PT"]))

algo = collect(tree["Jet/Jet.PT"])

# plot algo hist
histogram(algo, title="Jet PT", xlabel="PT", ylabel="Frequency")