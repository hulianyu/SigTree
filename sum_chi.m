function [pi,pval] = sum_chi(X,m,q)
pi = (X(:,m)~=q) + 1;
X_drop = X;
X_drop(:,m) = [];
M_drop = size(X_drop,2); %Attribute Number
dfs = 0;
chis = 0;
for m_drop=1:M_drop
    x = X_drop(:,m_drop);
    df = (length(unique(x))-1);
    if df~=0
        [tbl,chi,~] = crosstab(x,pi);
        %%
        % Check if the table is 2x2
        % if size(tbl, 1) == 2 && size(tbl, 2) == 2
        if size(tbl, 1) == 2
            % Calculate expected frequencies
            n = sum(tbl(:));  % Total sample size
            expected = zeros(size(tbl));
            for i = 1:2
                for j = 1:2
                    expected(i, j) = sum(tbl(i, :)) * sum(tbl(:, j)) / n;
                end
            end
            % Check if any expected frequency is less than 5
            if any(expected(:) < 5)
                % Apply Yates correction if any expected frequency is less than 5
                a = tbl(1, 1);
                b = tbl(1, 2);
                c = tbl(2, 1);
                d = tbl(2, 2);
                chi = (n * (abs(a * d - b * c) - n / 2)^2) / ((a + b) * (c + d) * (a + c) * (b + d));
            end
        end
        %%
        chis = chis + chi;
        dfs = dfs + df;
    end
end
pval = chi2cdf(chis,dfs,'upper');
end