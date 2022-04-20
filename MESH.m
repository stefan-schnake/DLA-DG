classdef MESH < handle
    %Stores mesh
    
    properties
        a;
        b;
        N;
        nodes;
        dx;
        jac;
    end
    
    methods
        function mesh = MESH(a,b,N)
            mesh.a = a;
            mesh.b = b;
            mesh.N = N;
            mesh.nodes = a:(b-a)/N:b;
            mesh.dx  = mesh.nodes(2)-mesh.nodes(1);
            mesh.jac = mesh.dx/2;
        end
    end
end

