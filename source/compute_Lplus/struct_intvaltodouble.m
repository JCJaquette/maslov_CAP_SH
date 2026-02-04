function [out,radii] = struct_intvaltodouble(struct_in)

    if isstruct(struct_in)

        out = struct();
        fields = fieldnames(struct_in);

        for k = 1:numel(fields)
            f = fields{k};
            [out.(f),radii.(f)] = struct_intvaltodouble(struct_in.(f));
        end

    else
        out = mid(struct_in);
        radii = rad(struct_in);
    end

end
