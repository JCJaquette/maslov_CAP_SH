function out = struct_doubletointval(struct_in)

    if isstruct(struct_in)

        out = struct();
        fields = fieldnames(struct_in);

        for k = 1:numel(fields)
            f = fields{k};
            out.(f) = struct_doubletointval(struct_in.(f));
        end

    else
        out = intval(1)*struct_in;
    end

end