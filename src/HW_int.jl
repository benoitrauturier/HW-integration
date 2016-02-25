
module HW_int


	# question 1 b)
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# here are some functions I defined for useage
	# in several sub questions

	# demand function

	# gauss-legendre adjustment factors for map change

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply,


	# weighted sum for integration from the slides.

	#Defining the function of interest, will be called in the subsequent questions
	function demand(p)
	  return 2.0.*p.^(-0.5)
	end

	function question_1b(n)
		#First we get the weights and points from the FastGaussQuadrature package
		#The variable name is PaW as "Points and Weights"
		PaW=gausslegendre(n)

		################COMPUTING THE SURPLUS####################################
		#Applying the formula over all the integration points. res is the variable which will contain in fine the change in consummer surplus.
		res=0
		for i in 1:n
			res+=demand(1.5*PaW[1][i]+2.5)*PaW[2][i]
		end
		#Now transforming the result to correspond to the given formula
		res*=1.5
		#Adding and substracting the rectangles
		res-=2
		##############PRINTING THE RESULT#######################################
		println("The estimated change in consummer surplus using Gauss-Legendre method is: $res")
		######################PLOT THE INTEGRATION POINTS########################
		figure(1)
		plot(1.5*PaW[1]+2.5,demand(1.5*PaW[1]+2.5),"bo")
		title("Gauss-Legendre")
		xlabel("Prices")
		ylabel("Quantities")
	end


	function question_1c(n)
		#setseed(20162302)
		#To generate random numbers between 1 and 4, one take random numbers from a uniform distribution between 0 and 1 and multiply it by 3 and add 1
		Points=rand(n)*3+1
		##################COMPUTING THE SURPLUS####################################
		#Compute the value of the function for all points
		res=0
		for i in 1:n
			res+=demand(Points[i])
		end
		#Average the sum of the values and multiply the result by the size of the interval
		res=3*res/n
		#Adding and substracting the squares
		res-=2
		##############PRINTING THE RESULT#######################################
		println("The estimated change in consummer surplus using Monte-Carlo integration method is: $res")
		######################PLOT THE INTEGRATION POINTS########################
		figure(2)
		plot(Points,demand(Points),"bo",color="red")
		title("Monte-Carlo")
		xlabel("Prices")
		ylabel("Quantities")
	end


	function question_1d(n)
		#Generating the sobol sequence; Not sure I completely understood that objec
		sob=SobolSeq(1,1,4) #This generate a sequence of points in [1,4]
		sobseq=ones(n)
		for i in 1:n
			sobseq[i]=next(sob)[1]
		end
		res=0
		for i in 1:n
			res+=demand(sobseq[i])
		end
		#Average the sum of the values and multiply the result by the size of the interval
		res=3*res/n
		#Adding and substracting the squares
		res-=2
		##############PRINTING THE RESULT#######################################
		println("The estimated change in consummer surplus using Pseudo Monte-Carlo integration method is: $res")
		######################PLOT THE INTEGRATION POINTS########################
		figure(3)
		plot(sobseq[:],demand(sobseq[:]),"bo",color="pink")
		title("Pseudo Monte-Carlo")
		xlabel("Prices")
		ylabel("Quantities")
	end


	function question_2a(n)
		#Intialize the dd function that will be needed in the code
		function dd(p,t1,t2)
        exp(t1)/p .+ exp(t2) .* p.^(-0.5) - 2
		end
############################GRID CONSTRUCTION##################################
		PaW=gausshermite(n)
		#Getting the Gauss Hermite Points and Weights
		#Constructing the weights vector
		Weights=kron(PaW[2],PaW[2])
		#Constructing the grid (n*n points in total)
		nodes=[]
		push!(nodes,repeat(PaW[1],inner=[n],outer=[1]))
		push!(nodes,repeat(PaW[1],inner=[1],outer=[n]))
		#Here we just change the type of the variable nodes to a simple array, instead of an array of array which would yield erros in the subsequent parts of the code
		nodes2=[nodes[1] nodes[2]]
##########################GRID ADJUSTMENT######################################
		Sigma=[0.02 0.01;0.01 0.01]
		Omega = chol(Sigma,Val{:L})
		adj_nodes=Omega*transpose(nodes2)+zeros(2,n*n)
###############CALCULATION OF THE MEAN AND THE VARIANCE########################
		meanp=0
		meanp2=0
#		'''
#		Here the principle is that at each loop, we define a new "price function",
#		which is not closed form. Solving this thanks to the fzero function, we get
#		p for a particular point on the grid. Then we weight it. Finally, we add it
#		to the previous "weighted prices". We do the same for p^2.
#		'''
		for i in 1:n*n
			function price(p)
				dd(p,adj_nodes[1,i],adj_nodes[2,i])
			end
			#Note that this method yield approximations of the root. A few time the function price(p) returns a value different from zero. All the time the error seems to be neglibible (lower than 0.01)
			p=fzero(price,[0,10000000])
			#Here we assume that there is something unclear on slide #39, i.e. we devide by pi
			meanp2+=(p.^2).*Weights[i]./pi
			meanp+=p.*Weights[i]./pi
		end
		#Compute the variance with the formula VAR[X]=E[X^2]-E[X]^2
		Varp=meanp2-meanp^2
#########################PRINTING RESULTS###############################################
		print("The estimated expectation of the price with Gauss-Hermite grid is: $meanp \n")
		print("The estimated variance of the price with Gauss-Hermite grid is: $Varp \n")
##########################PLOTTING THE GRID#####################################
		figure(4)
		plot(adj_nodes[1,:],adj_nodes[2,:],"bo")
		title("Gauss-Hermite")
		xlabel("Theta 1")
		ylabel("Theta 2")
	end


	function question_2b(n)
		function dd(p,t1,t2)
				exp(t1)/p .+ exp(t2) .* p^-0.5 - 2
		end
		n=n*n #This is to have the same number of grid points with runall
		Sigma=[0.02 0.01;0.01 0.01]
		Law=MvNormal(zeros(2),Sigma)
		grid=rand(Law,n)
		meanp=0
		meanp2=0
#		'''
#		Here the principle is that at each loop, we define a new "price function",
#		which is not closed form. Solving this thanks to the fzero function, we get
#		p for a particular point on the grid. Then we weight it. Finally, we add it
#		to the previous "weighted prices". We do the same for p^2.
#		'''
		for i in 1:n
			function price(p)
				dd(p,grid[1,i],grid[2,i])
			end
			#Note that this method yield approximations of the roots. A few time the function price(p) returns a value different from zero. All the time the error seems to be neglibible (lower than 0.01)
			p=fzero(price,[0,1000000])
			meanp2+=(p[1].^2)./n
			meanp+=p[1]./n
		end
		#Compute the variance with the formula VAR[X]=E[X^2]-E[X]^2
		Varp=meanp2-meanp^2
#########################PRINTING RESULTS###############################################
		print("The estimated expectation of the price with Monte-Carlo grid is: $meanp \n")
		print("The estimated variance of the price with Monte-Carlo grid is: $Varp \n")
		######################PLOT THE GRID########################################
		figure(5)
		plot(grid[1,:],grid[2,:],"bo",color="green")
		title("Monte-Carlo")
		xlabel("Theta 1")
		ylabel("Theta 2")
	end


	# function to run all questions
	function runall(n=10)
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n)	# make sure your function prints some kind of result!
		question_1c(n)
		question_1d(n)
		println("")
		println("results of question 2:")
		q2 = question_2a(n)
		println(q2)
		q2b = question_2b(n)
		println(q2b)
		println("end of HW-integration")
	end

end
