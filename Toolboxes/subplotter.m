function subplotter(x,data)

 nSeries = size( data, 2 ) ;

 % - Build figure.
 figure() ;  clf ;
 set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
 % - Compute #rows/cols, dimensions, and positions of lower-left corners.
 nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
 rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
 colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
 rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
 % - Build subplots axes and plot data.
 for dId = 1 : nSeries
    rowId = ceil( dId / nCol ) ;
    colId = dId - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    plot( x(dId,:), data(:,dId), 'b' ) ;
    grid on ;
    xlabel( '\theta(t) [rad]' ) ;  ylabel( 'Anomaly [m]' ) ;
    title( sprintf( 'Time series %d', dId )) ;    
 end
 % - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
 text( 0.5, 0, 'My Nice Title', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end