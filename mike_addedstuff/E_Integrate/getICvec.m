function ICvec = getICvec(pulseIC,basisVec1,basisVec2)

    if dot(pulseIC,basisVec1) > 10^-3
        vec = basisVec1;
    else
        vec = basisVec2;
    end

    projection = dot(pulseIC,vec)/dot(vec,vec) * vec;

    ICvec = vec - projection;

end

