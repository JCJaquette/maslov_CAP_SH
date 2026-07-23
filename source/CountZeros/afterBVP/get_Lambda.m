function out = get_Lambda(params,str)

jacparams.nu = mid(params.nu);
jacparams.mu = mid(params.mu);
jacparams.scale = params.scale;

[~, values]= getJacEigs(0, jacparams);

out = values.(str);

end