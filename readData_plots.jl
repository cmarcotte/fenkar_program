using Plots, Printf, ProgressBars, Measures
gr()

x = rand(Float32,200,200,50,3);

zs = [1,10,20,30,40,50]
dx = 0.02

cb = []
push!(cb, scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), c=:Oranges, colorbar_title="\$u\$", colorbar_ticks=[0,1],
						xlims=(1,1), label="", framestyle=:none, margin=0mm, top_margin=0mm, bottom_margin=-10mm))
push!(cb, scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), c=:Purples, colorbar_title="\$v\$", colorbar_ticks=[0,1],
						xlims=(1,1), label="", framestyle=:none, margin=0mm, top_margin=10mm, bottom_margin=0mm))
push!(cb, scatter([0,0], [0,1], zcolor=[0,3], clims=(0,1), c=:Greens, colorbar_title="\$w\$", colorbar_ticks=[0,1],# bottom_margin=10mm,
						xlims=(1,1), label="", framestyle=:none, margin=0mm, top_margin=-1mm, bottom_margin=10mm))
cp = plot(cb...; layout=(3,1), margin=0mm, left_margin=0mm)#left_margin=-15mm, right_margin=0mm)#, bottom_margin=5mm, top_margin=5mm)

p = [plot() for j in 1:length(zs), i in 1:3]

traces = []
anim = Animation()
for t in ProgressBar(0:100)
	try
		fname = @sprintf("./out/fenkar_200_200_ 50_%4d.dat",t);
		#print("Reading from $(fname):\n")
		open(fname,"r") do io
		       read!(io, x);
		end
		for (o,z) in enumerate(zs)
			p[o,1] = heatmap(x[:,:,z,1]; 	c=:Oranges, clims=(0,1), #title="\$ z=$(round(z*dx,digits=3))\$ [cm]", top_margin=2mm,
									aspect_ratio=1.0, lims=[1,200], colorbar=false, xticks=:none, yticks=:none, margin=0mm)
			p[o,2] = heatmap(x[:,:,z,2];	c=:Purples, clims=(0,1),
									aspect_ratio=1.0, lims=[1,200], colorbar=false, xticks=:none, yticks=:none, margin=0mm)
			p[o,3] = heatmap(x[:,:,z,3]; 	c=:Greens, clims=(0,1), xlabel="\$ z=$(round(z*dx,digits=3))\$ [cm]", bottom_margin=10mm,
									aspect_ratio=1.0, lims=[1,200], colorbar=false, xticks=:none, yticks=:none, margin=0mm)
		end
		#hp = plot(p...; layout=(3,6), link=:all)
		hp = plot(p...; layout=(3,6), link=:all, margin=0mm)
		l = @layout [a b{0.07w}]
		newframe = plot(hp, cp; size=(1.1*200*length(zs), 200*3*1.1), margin=0mm, layout=l, plot_title="\$t=$(t)\$ [ms]")
		frame(anim, newframe)
		push!(traces, x[100,100,25,:])
	catch
	end
end
mp4(anim, "./heatmapsPlots.mp4", fps=10)


plot(2.0.*(0:100), transpose(reduce(hcat, traces)), xlabel="\$ t \$ [ms]")
savefig("./tracesPlots.pdf")
