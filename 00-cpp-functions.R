library(Rcpp)

cppFunction('NumericVector czf(NumericVector x) {
            int n = x.size();
            NumericVector out(n);
            for(int i = 0; i < n; ++i) {
            if(i > 0 && out[i - 1] == 0 && NumericVector::is_na(x[i])){
            out[i] = 0;
            } else {
            out[i] = x[i];
            }
            }
            return out;
            }')

cppFunction('NumericVector mark_bounds(NumericVector x, int steps, NumericVector lower_bound){
            int n = x.size();
            NumericVector out(n);
            for(int i = 0; i < n; ++i) {
            if(i > 0 && NumericVector::is_na(x[i]) && !NumericVector::is_na(x[i - 1]) && x[i - 1] <= lower_bound[0] ){
            bool future_all_NA = true;
            for(int j=i+1; j <= i + steps; j++){
            future_all_NA = future_all_NA && NumericVector::is_na(x[j]);
            }
            if(future_all_NA){
            out[i] = 0;
            } else {
            out[i] = x[i];
            }
            } else {
            out[i] = x[i];
            }
            }
            return out;
            }')

cppFunction('NumericVector locf(NumericVector x, int steps){
            int n = x.size();
            NumericVector out(n);
            for(int i = 0; i < n; ++i) {
            if(i > 0 && NumericVector::is_na(x[i]) && !NumericVector::is_na(x[i - 1])){
            int start = i;
            for(; i < start + steps; i++){
            if(i< n && NumericVector::is_na(x[i])){
            out[i] = x[start - 1];
            } else {
            break;
            }
            }
            i--;
            } else {
            out[i] = x[i];
            }
            }
            return out;
            }')

cppFunction('CharacterVector locfStr(CharacterVector x, int steps){
            int n = x.size();
            CharacterVector out(n);
            for(int i = 0; i < n; ++i) {
            if(i > 0 && CharacterVector::is_na(x[i]) && !CharacterVector::is_na(x[i - 1])){
            int start = i;
            for(; i < start + steps; i++){
            if(i< n && CharacterVector::is_na(x[i])){
            out[i] = x[start - 1];
            } else {
            break;
            }
            }
            i--;
            } else {
            out[i] = x[i];
            }
            }
            return out;
            }')

cppFunction('DataFrame modifySingleColumn(DataFrame df){
            DataFrame out(df);
            IntegerVector v = out["speed"];
            for ( int i =0 ; i < v.size(); i++){
            v[i]--;
            }
            out["speed"] = v;
            return out;
            }')

rapid_imputation <- function(x, steps, lower = -Inf) {
  result <- mark_bounds(x, steps = steps, lower = lower)
  result <- locf(result, steps = steps)
  result <- czf(result)
  return(result)
}