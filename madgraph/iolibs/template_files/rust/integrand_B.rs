/*
 * This is the rust implementation of the integrand %(integrand_name)s
 *
 */

## if (my_conditional_variable > 3) {
    println!("A");
## } else {
    println!("B");
## if (my_nested_conditional_boolean) {
    println!("B1");
## } else {
    println!("B2");
## }
## }
    
    println!("This is inconditional");

## if (any(var=='PrintTHIS' for var in myContextList)) {
    println!("Some statement dynamically assigned:")
    %(a_dynamical_statement)s
## }
