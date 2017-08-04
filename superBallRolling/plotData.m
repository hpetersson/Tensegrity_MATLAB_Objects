function plotData(a1, a2)
    global stringTensionDatastore;
    global stringTensionCmdDataStore;
    global restLenDataStore;
    
    plot(stringTensionDatastore); hold on;
    plot(stringTensionCmdDataStore, 'k');
    title('Tensions');
    hold off;
    figure();
    plot(restLenDataStore);
    title('Rest lengths');
    
end

