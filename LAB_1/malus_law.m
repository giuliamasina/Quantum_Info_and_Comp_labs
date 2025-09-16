function result = malus_law(I0, theta, theta0, offset)
    result = I0 * cos(theta - theta0).^2 + offset;
end
