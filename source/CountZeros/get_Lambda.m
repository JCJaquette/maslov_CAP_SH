function out = get_Lambda(params)

jacparams.nu = mid(params.nu);
jacparams.mu = mid(params.mu);
jacparams.scale = params.scale;

[~, values]= getJacEigs(0, jacparams);
out = values.u;

end