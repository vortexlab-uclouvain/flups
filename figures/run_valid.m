function run_valid(valid_case,order)

    err = load( ['../data/validation_' valid_case '.err']);
    N = err(:,1);
    err2 = err(:,2);
    erri = err(:,3);

    figure(); hold all;
    plot(N,err2,'.-');
    plot(N,erri,'.-');
    set(gca,'Xscale','log','Yscale','log');
    grid on;
    box on;
    
    plot(N,N.^(-order)*erri(end)/N(end)^(-order),'--k');
    plot(N,N.^(-order)*err2(end)/N(end)^(-order),'--k');
    
    set(gca,'XTick',N);
    title(valid_case);
    legend('err2','erri');
    
 
end
