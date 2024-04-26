function plot_edges(connspring,x,ee)
      hold on; 
      for e = 1:size(  connspring,1) % loop over all string elements
    
        % extract node position from current string element   
        iie  =   connspring(e, :); % nodes of the curren  t element e
        xe   =   x(iie);      % coordinates of these nodes
        plot(dot(xe,ee(1)),dot(xe,ee(2)),"LineWidth",1.5)
      end
      hold off;
end