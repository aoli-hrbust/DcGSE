<<<<<<< HEAD
function ACC = Accuracy(C,gt)
 C = bestMap(gt,C);
 ACC = length(find(gt == C))/length(gt);
=======
function ACC = Accuracy(C,gt)
 C = bestMap(gt,C);
 ACC = length(find(gt == C))/length(gt);
>>>>>>> e1248d18d20a74f08f96760f80be9bf1df35eab1
end