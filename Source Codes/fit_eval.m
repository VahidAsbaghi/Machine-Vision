function [s_prim]=fit_eval(out_taw,s,ci,ci1)
if ci<=ci1
    s_prim=s^out_taw;
else
    s_prim=1-(1-s)^out_taw;
end
end
