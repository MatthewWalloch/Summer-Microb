# roulette wheel selection 
function fortune_wheel(weights::Vector{Float64})

  accumulation = cumsum(weights)
  p = rand() * accumulation[end]
  index = 0
  while index<length(accumulation)
    index+=1
    if accumulation[index] > p
      print(index)
      return index
    end
  end

end
