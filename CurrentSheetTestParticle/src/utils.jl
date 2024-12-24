function distance(A, B)
    n, m = size(A, 1), size(B, 1)
    result = zeros(n)
    Threads.@threads for i = 1:n
        minSqDist = Inf
        @inbounds for j = 1:m
            dx = A[i, 1] - B[j, 1]
            dy = A[i, 2] - B[j, 2]
            dz = A[i, 3] - B[j, 3]
            sqDist = dx * dx + dy * dy + dz * dz
            if sqDist < minSqDist
                minSqDist = sqDist
            end
        end
        result[i] = minSqDist
    end
    return minimum(result)
end

distance(sol1::ODESolution, sol2::ODESolution) = distance(sol1[1:3, :]', sol2[1:3, :]')

function get_gc_func(B)
    param = prepare(E, B, species=User)
    get_gc(param)
end