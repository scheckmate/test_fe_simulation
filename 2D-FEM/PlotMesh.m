function [ ] = PlotMesh( Nodes, Elements, c, L, Color, ax)
%Plot the mesh given by "Nodes" and "Elements".
%"c" and "L" are used to get some space between the mesh and the boundary of
%the figure.
%If "Color" is zero the elements are colored white. Else "Color" should be a
%vector of size number of elements and the elements are colored according
%to the values of "Color".
if(Color == 0)
    patch('Vertices',Nodes,'Faces',Elements,'FaceColor','w','Parent',ax);
else
    patch('Vertices',Nodes,'Faces',Elements,'FaceColor','flat','FaceVertexCData',Color,'Parent',ax);
    colorbar;
    %h = colorbar;
    %set(h,'YTickLabel',{num2str(min(Color)), num2str(max(Color))});
end
axis equal;
axis([-1 L+5 -c-5 c+1]);
xlabel('X');
ylabel('Y');

end

