function out = flattenstruct(s, prefix)

    out = struct();
    fields = fieldnames(s);

    for i = 1:numel(fields)
        f = fields{i};
        val = s.(f);

        if isempty(prefix)
            name = f;
        else
            name = [prefix '_' f];
        end

        if isstruct(val)
            sub = flatten_struct(val, name);
            subfields = fieldnames(sub);
            for j = 1:numel(subfields)
                out.(subfields{j}) = sub.(subfields{j});
            end
        else
            out.(name) = val;
        end
    end

end