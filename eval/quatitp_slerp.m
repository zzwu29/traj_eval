function Mo = quatitp_slerp(Mi, ti, to, idx)

s   = (to - ti(idx(:, 1)))./(ti(idx(:, 2)) - ti(idx(:, 1)));

% dM = zeros(size(idx,1), size(Mi,1));
% for i=1:size(Mi,1)
%     dM=quatmultiply(quatinv(Mi(idx(:, 1), :)), Mi(idx(:, 2), :));
% end

dM = quatmultiply(quatinv(Mi(idx(:, 1), :)), Mi(idx(:, 2), :));
% dM  = Mi(idx(:, 2), :) - Mi(idx(:, 1), :);
sdM = [dM(:, 1).^s, dM(:, 2).^s, dM(:, 3).^s, dM(:, 4).^s];

Mo  = quatmultiply(Mi(idx(:, 1), :), sdM);

end