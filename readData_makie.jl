using GLMakie, Printf, ProgressBars

	times = 0:500
	
	x = rand(Float32,200,200,50,3);
	
	U = Observable(x[:,:,:,1])
	V = Observable(x[:,:,:,2])
	W = Observable(x[:,:,:,3])

	fig = Figure(resolution = (1600, 600))
	axs = [Axis3(fig[1,i], aspect=:data, elevation=0.75, azimuth=2.65, perspectiveness=0.2) for i in 1:3]
	volume!(axs[1],U,algorithm=:mip,colormap="Oranges",alpha=0.05,transparence=true)
	volume!(axs[2],V,algorithm=:mip,colormap="Purples",alpha=0.05,transparence=true)
	volume!(axs[3],W,algorithm=:mip,colormap="Greens",alpha=0.05,transparence=true)
	for ax in axs
		hidedecorations!(ax)
	end
	Colorbar(fig[2,1]; vertical = false, flipaxis = false, colorrange=(0,1), colormap="Oranges", label=L"$u$")
	Colorbar(fig[2,2]; vertical = false, flipaxis = false, colorrange=(0,1), colormap="Purples", label=L"$v$")
	Colorbar(fig[2,3]; vertical = false, flipaxis = false, colorrange=(0,1), colormap="Greens",  label=L"$w$")
	#tit = Label(fig[0, :], text = L"$t = 0$ [ms]", textsize = 32)
	resize_to_layout!(fig)

	record(fig, "./volumeMakie.mp4", ProgressBar(times); framerate=20) do t
		try
			fname = @sprintf("./out/fenkar_200_200_ 50_%4d.dat",t)
			#print("Reading from $(fname):\n")
			open(fname,"r") do io
			       read!(io, x);
			end
			#tit.text[] = L"$t = %$(2.0*t)$ [ms]"
			U[] = x[:,:,:,1]
			V[] = x[:,:,:,2]
			W[] = x[:,:,:,3]
		catch
		end
	end



	Zs = [1,10,20,30,40,50]
	us = [Observable(x[:,:,z,1]) for z in Zs]
	vs = [Observable(x[:,:,z,2]) for z in Zs]
	ws = [Observable(x[:,:,z,3]) for z in Zs]

	fig = Figure(resolution = (1450, 750))
	axs = [ Axis(fig[i, j], width = 200, height = 200) for i in 1:3, j in 1:length(Zs) ]

	for ax in axs
		hidedecorations!(ax)
	end
	
	for j in 1:length(Zs)
		heatmap!(axs[1,j],us[j],colormap="Oranges",colorrange=(0,1))
		heatmap!(axs[2,j],vs[j],colormap="Purples",colorrange=(0,1))
		heatmap!(axs[3,j],ws[j],colormap="Greens", colorrange=(0,1))
		axs[1,j].title = "$(round(Zs[j]*0.02; digits=3)) [cm]"	# makie handles latex very confusingly; only at axis creation time, only subset of latex, not using latex explicitly
		#axs[1,j].xlabel = L" "
		#axs[2,j].xlabel = L" "
		#axs[3,j].xlabel = L"$z = %$(round(Zs[j]*0.02; digits=1))$ [cm]"
	end
	Colorbar(fig[1,length(Zs)+1]; colorrange=(0,1), colormap="Oranges", label=L"u")
	Colorbar(fig[2,length(Zs)+1]; colorrange=(0,1), colormap="Purples", label=L"v")
	Colorbar(fig[3,length(Zs)+1]; colorrange=(0,1), colormap="Greens",  label=L"w")
	#tit = Label(fig[0,:], text = L"$t = 0$ [ms]", textsize = 32)
	resize_to_layout!(fig)

	traces = []
	# create video directly; blurrier, faster (if using GLMakie backend)
	record(fig, "heatmapsMakie.mp4", ProgressBar(times); framerate=10) do t
		fname = @sprintf("./out/fenkar_200_200_ 50_%4d.dat",t)
		#print("Reading from $(fname):\n")
		open(fname,"r") do io
		       read!(io, x);
		end
		#tit.text[] = L"$t = %$(2.0*t)$ [ms]"
		for (j,z) in enumerate(Zs)
			us[j][] = x[:,:,z,1]
			vs[j][] = x[:,:,z,2]
			ws[j][] = x[:,:,z,3]
		end
		push!(traces, x[100,100,25,:])
	end
	traces = reduce(hcat,traces)

	var = ["u","v","w"]
	fig = Figure(resolution = (1200,400))
	ax = Axis(fig[1,1], width=1100, height=300; xlabel=L"$ t $ [ms]")
	for n in 1:3
		lines!(ax, 2.0.*times, traces[n,:], label=L"%$(var[n])")
	end
	fig[0,1] = Legend(fig, ax, "", framevisible = false, nbanks=3)
	resize_to_layout!(fig)
	save("./tracesMakie.png", fig)

