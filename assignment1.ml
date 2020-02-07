(*
		Author: Daman Arora
		Date: February 1, 2020
		Course: COL226(Programming Languages)
		Purpose:
			This code is built to perform certain operation on matrices and vectors
			This contains a working implementation of Gauss-Jordan elimination.
		Caution:
			This code works on floats which won't be suitable for high-precision 
			computation for implementing a Gauss-Jordan elimination algorithm. This
			has been designed keeping in mind speed and efficiency rather than precision.
			For precision, use a format with higher precision than 'float'
*)

type vector = float list;;
type matrix = float list list;;

exception InvalidInput;;
exception UnequalVectorSize;;
exception UnequalMatrixShape;;
exception IncompatibleMatrixShape;;
exception SingularMatrix;;

(*		
	vdim works recursively on the elements of the list. 
	Base case -> empty list []
	Complexity Analysis-> 
		By the recursive defn of the function, we have 
			T(n) = T(n - 1) + 1
		Hence, the complexity is O(n) where n is the length of the list.
	Test cases:
		1) []
		2) [3.; 4.; 9.4; 6.]
		3) [3.]
*)
let rec vdim (v:vector): int =
	match v with
  		| [] -> 0
  		| x :: xs -> 1 + vdim(xs);;

(*
	mkzerov appends a 0. to a list and goes on iteratively.
	Complexity Analysis->
		This also follows the same logic of vdim and hence is O(n)
	Test cases:
		1) 0
		2) 5
*)

let rec mkzerov (n:int): vector =
	if n = 0 then 
		[]	
  	else 
  		0.0 :: (mkzerov(n - 1) : vector);;

(*
	iszerov also works on the same principle. It is O(n)
	Possible improvement:
		This algorithm does better on an amortized level is the 
		language compiler implements boolean short-circuiting.
	Test cases:
		1) []
		2) [0.]
		3) [0.; 0.; 1.]
*)

let rec iszerov (v:vector): bool = match v with 
	| [] -> true
	| x :: xs -> (x = 0.0) && iszerov(xs);;

(*
	The algorithm has a complexity of O(min(l1, l2)) where
	l1 and l2 are the dimensions of the two vectors.
	Test case:
		1) [] []
		2) [] [4.; 5.]
		3) [1.2; 9.] [2.4; 5.]
*)

let rec addv (v1:vector) (v2:vector): vector = match v1, v2 with
	| [], [] -> []
	| [], y :: ys -> raise UnequalVectorSize;
	| x :: xs, [] -> raise UnequalVectorSize;
	| x :: xs, y :: ys -> (x +. y) :: addv xs ys;;

(*
	scalarmultv is O(n) where n is the dimension of the vector
	Test case:
		1) 5. []
		2) 8. [0.; 5.; 6.]
		3) 0. [5.6]
*)

let rec scalarmultv (c:float) (v:vector): vector = match v with
	| [] -> []
	| x :: xs -> (x *. c) :: (scalarmultv c xs : vector);;


(*
	dotprodv is O(min(l1, l2)) where l1 and l2 are the dimensions
	of the two vectors
	Test case:
		1) [] []
		2) [] [5.]
		3) [4.; 5.; 6.] [-1.; 3.; 0.]
*)
let rec dotprodv (v1:vector) (v2:vector): float = match v1, v2 with
	| [], [] -> 0.0
	| [], y :: ys -> raise UnequalVectorSize;
	| x :: xs, [] -> raise UnequalVectorSize;
	| x :: xs, y :: ys -> (x *. y) +. dotprodv xs ys;;

(*
	mdim also performs a check on whether all rows have equal lengths or not using
	the 'isValid' function.
	Complexity Analysis ->
		isValid is O(num) where num is the number of elements in the entire matrix.
		For a valid matrix is would be O(m * n)

		mdim is also O(m * n)
	Test case:
		1) []
		2) [[]]
		3) [[2.3; 4.5]; [4.6; 7.1]; [8.1; 9.0]]
		4) [[1.]; [2.]; [3.; 4.]]
*)

let rec isValid (mat : matrix) dim = match mat with
	| [] -> true
	| x :: xs -> 
		if(List.length x = dim) then 
			isValid xs dim
		else
			false;;

let rec mdim (m:matrix): int*int = match m with 
	| [] -> (0, 0)
	| x :: xs ->
		let colLength = vdim x in 

		if(isValid m (List.length x)) then 
			(List.length m, colLength)
		else	
			raise UnequalMatrixShape;;

(*	mkzerom is O(m * n). It creates m lists of n zeros using the mkzerov subroutine
	Test cases:
		1) 0 0
		2) 1 0
		3) 0 1
		4) 3 4
*)


let rec mkzerom (m_:int) (n_:int): matrix =	
	if m_ = 0 then 
		[]
	else 
		mkzerov(n_) :: mkzerom (m_ - 1) n_;;

(*	iszerom is O(m * n). It checks the shape-correctness of the matrix as well.
	Test cases:
		1) []
		2) [[0.; 0.]; [0.; 0.]]
		3) [[0.; 0.]; [1.]]
*)

let rec isvalidzero (mat : matrix) dim= match mat with
	| [] -> true
	| x :: xs -> 
		if (List.length x) = dim then 	
			if (iszerov x) then 
				isvalidzero xs dim
			else 
				false
		else
			raise UnequalMatrixShape;;

let rec iszerom (m:matrix): bool = match m with
	| [] -> true
	| x :: xs -> isvalidzero m (List.length x);;


(*
	mkunitm creates matrices iteratively using the mkunitIter function
	It uses the mkunitv subroutine which creates a vector of the form
	[...0.; 0.; .... 1.; 0.; 0....]
	Complexity is O(n * n)
	Test cases:
		1) -1
		2) 0
		3) 10
*)
let rec mkunitv i n size = 	
	if(i = size) then 
		[]
	else if(i = n) then
		1.0 :: (mkunitv (i + 1) n size : vector)
	else 
		0.0 :: mkunitv (i + 1) n size;;
							
let rec mkunitIter i m = 
	if(i = m) then 
		[]
	else
		(mkunitv 0 i m) :: mkunitIter (i + 1) m;;

let rec mkunitm (m_:int): matrix =
	if(m_ >= 0) then 
		mkunitIter 0 m_
	else 
		raise InvalidInput;;


(*
	isunitm is O(n * n). It checks whether each row is correct or not. 
	It assumes dimension-correctness. 
	Test cases:
		1) [[]]
		2) mkunitm 10
		3) [[1.2]]
*)
let rec isonev v i n=   
	match v with
		| [] -> true
		| (x :: xs) ->  
			if(i = n) then 
				(x = 1.0) && isonev xs (i + 1) n
			else 
				(x = 0.0) && isonev xs (i + 1) n;;

let rec checkunit i (mat : matrix) dim = 
	match mat with
		| [] -> 
			true
		| x :: xs -> (isonev x 0 i) && checkunit (i + 1) xs dim;;

let rec isunitm (m:matrix): bool =
	checkunit 0 m (List.length m);;


(*
	addm is O(n * n). It adds rows using the addRows function. 
	It checks whether each row is correct or not. 
	Test cases:
		1) [[1.]] [[1.2]]
		2) [[1.2]] [[1.2; 5.6; 8.]]
		3) mkunit 8 mkunit 8
*)

let rec addRows (v1:vector) (v2:vector): vector = match v1, v2 with
	| [], [] -> []
	| [], y :: ys -> raise UnequalMatrixShape;
	| x :: xs, [] -> raise UnequalMatrixShape;
	| x :: xs, y :: ys -> (x +. y) :: addRows xs ys;;

let rec addm (m1:matrix) (m2:matrix): matrix =
	match m1, m2 with
		| [], [] -> []
		| [], xs -> raise UnequalMatrixShape
		| xs, [] -> raise UnequalMatrixShape
		| x :: xs, y :: ys -> addRows x y :: addm xs ys;;


(*	scalarmultm has a complexity of O(m * n) where m and n are 	
	the dimensions of the matrix. It uses the scalarmultv routine.
	Test cases:
		1) 2. [[3.4; 5.6; 7.]; [0.; 0.; 0.]]
		2) 0. [[]]
*)
let rec scalarmultm (c:float) (m:matrix): matrix =  
	match m with
		| [] -> []
		| x :: xs -> scalarmultv c x:: scalarmultm c xs;;

(* 
	append is a function that elements of v to the end of rows
	of m matrix. 
	Eg -> If we append (5, 6) to 
		[1, 2]
		[3, 4], we get
		[1, 2, 5]
		[3, 4, 6]

	It has complexity O(r * c) where r and c are the number of rows and 
	columns in the matrix m;


	transm uses this to invert each row in each iteration.
	Complexity analysis: 
		For the ith row, number of steps requires is O(i * c) where c is the
		number of columns. For doing this for r rows, complexity achieved is 
		O(r ^ 2 * c).

	Improvement:
		This could be done in O(r * c) if we dealt with arrays instead of lists.

	Test cases:
		1) mkunit 10
		2) mkzerom 10 11
*)



let rec append (m : matrix) (v : vector)= 
	match m, v with
		| [],  [] -> []
		| [], y :: ys -> [y] :: append [] ys
		| x :: xs, [] -> []
		| x :: xs, y :: ys -> (x @ [y]) :: append xs ys;;

let rec transIter (m : matrix) (t : matrix) = 
	match m with
		| [] -> t
		| x :: xs -> transIter xs (append t x);;

let rec transm (m:matrix): matrix =  transIter m [];;



(* 	multm performs a row-by-row multiplication using the transpose of the 
	second matrix. This is computationally more feasible since row traversal
	is less expensive. 
	Complexity Analysis:	
		multRowByRow takes O(r1 * r2 * c * c) computations where r1 and c are the dimensions
		of the first matrix and c and r2 are the dimensions of the second matrix;

		multm only performs dimension checking and transpose besides multRowByRow. 
		Dimension checking has complexity O(c * (r1 + r2)).
		Transpose takes time O(r2 ^ 2 * c).
		Hence multm can be done in
			O(r2^2*c + r1*c + r2*c + r1*r2*c^2)

	Test cases:
		1) (mkzerom 2 3) (mkzero 3 5)
		2) (mkzerom 2 3) (mkzero 4 5)
		3) [[1.; 2.]; [3.; 4.]] [[5.; 6.]; [7.; 8.]]
*)

let rec multRowByRow (mat1 : matrix) (mat2 : matrix) = 
	match mat1 with
		| [] -> []
		| x :: xs -> List.map (dotprodv x) mat2 :: (multRowByRow xs mat2);;

let rec multm (m1:matrix) (m2:matrix): matrix = 
	let (_, c1) = mdim m1 in 
	let (c2, _) = mdim m2 in 
	if(c1 = c2) then
		multRowByRow m1 (transm m2)
	else
		raise IncompatibleMatrixShape;;


(*
	For calculating inverse, we first convert the matrix into upper-triangular form.
	After that, we multiply the elements of the diagonal and them compute the determinant.
	For this, Gauss-Jordan elimination is used.
	Firstly row-operations are done to ensure upper-triangular form.



	Algorithm:	
		Let the matrix be m with size s
		
		for each col c: 	
			find elem e(pivot) such that |e| is maximum in m[c, c]...m[c, s]
			if(e = 0)
				proceed to the next column
			else
				swap pivotRow and m[c]. (maintain count of swaps)
				Scale all rows below m[c] such that m[[i][c] = pivot or 0.
				Subtract m[c] from each column
				Store product for each scaling factor. 
		Compute product of diagonal entries
		Divide by product of scaling factors.


	Purpose of helper functions:

		findMaxIndexInColumn: This computes the maximum value index in a column of a matrix.
		findPivotIndex: Trims all the part above the current row and returns pivot index.
		swap, swapRows: swaps two rows r1 and r2.
		subtract: Computes difference of two vectors.
		scaleFactor, mult: Determine scale factors during scaling process.
		scaleRowsBelow: Multiplies all rows with corresponding scale factors.
		generateZerosBelow: Subtracts pivot row from each row.
		generateGaussForm: Generates the Gaussian form of the matrix.
		diagProd: Computes the product of diagonal entries.

	Complexity Analysis:	
		findMaxIndexInColumn:  
			O(r * col) where r is no of rows in the matrix 
			and col is the column
		findPivotIndex:
			O(col ^ 2) where col is the column number
		swap, swapRows: 
			O(r) where r is the number of rows
		subtract: 
			O(n) where n is the dimension of the vector
		scaleRowsBelow:
			O(n^2 - col^2) where n is the dimension of the matrix and col 
			is the column number
		generateZerosBelow:
			Identical to scaleRowsBelow O(n^2 - col^2)
		generateGaussForm:
			Scaling for column i takes time O(n ^ 2) as a rough estimate.
			Hence for each column, we require O(n^2) time.
			->For the entire matrix we require O(n^3) time.
			Computing diagonal product takes O(n^2).
			
		Hence, complexity is O(n^3).
		Test cases:
			1) [[0.]]
			2) [[1.; 2.]; [2.; 4.]]
			3) [[1.; 3.; 0.]; [2.; 5.; 1.]; [1.; 2.; 1.]]		
			4) [[1.; 3.; 2.; 0.]; [0.; -3.; -4.; 2.]; [-1.; -2.; 0.; 4.]; [-4.; 3.; 1.; 4.]]					
			5) scalarmultm 4. (mkunitm 100)
*)

let absf(f : float) = if f > 0.0 then f else 0.0 -. f;;

let rec findMaxIndexInColumn (mat:matrix) col currentIndex prevIndex prevValue = match mat with
	| [] -> prevIndex
	| x :: xs ->	
		let curr = absf(List.nth x col) in
		if(curr > prevValue) then 
			findMaxIndexInColumn xs col (currentIndex + 1) currentIndex curr
		else
			findMaxIndexInColumn xs col (currentIndex + 1) prevIndex prevValue;;

let rec findPivotIndex (mat : matrix) idx col = 
	match mat with
		| [] -> raise UnequalMatrixShape		
		| x :: xs -> 	
			if(idx = col) then
				findMaxIndexInColumn mat col col col (absf(List.nth x col))
			else
				findPivotIndex xs (idx + 1) col;;

let rec swap (mat:matrix) r1 r2 row1 row2 = 
	match mat with
		| [] -> []
		| x :: xs -> 	
			if(r1 = 0) then 
				row2 :: swap xs (r1 - 1) (r2 - 1) row1 row2
			else if(r2 = 0) then 
				row1 :: swap xs (r1 - 1) (r2 - 1) row1 row2
			else
				x :: swap xs (r1 - 1) (r2 - 1) row1 row2;;

let rec swapRows (mat:matrix) r1 r2 = 	
	let row1 = List.nth mat r1 in
	let row2 = List.nth mat r2 in
	if(r1 = r2) then 
		(mat, 1.0)
	else
		(swap mat r1 r2 row1 row2, (-1.0));;

let rec subtract (v1:vector) (v2:vector) = 
	match v1, v2 with
		| [], [] -> []
		| [], y :: ys -> raise UnequalMatrixShape
		| x :: xs, [] -> raise UnequalMatrixShape
		| x :: xs, y :: ys -> (x -.y) :: (subtract xs ys);;

let scaleFactor a b = if(b = 0.0) then 1.0 else a /. b;;

let mult f a = f *. a ;;

let rec scaleRowsBelow (mat:matrix) idx col pivot =
	match mat with
		| [] -> ([], 1.0)
		| x :: xs -> 	
			if(pivot = 0.0) then 
				(mat, 1.0)
			else
				let denom = List.nth x col in
				let sf = scaleFactor pivot denom in
				let multiplier = mult sf in
				let (next, prodNext) = scaleRowsBelow xs (idx + 1) col pivot in
				if(idx < col) then
					(x :: next, prodNext)
				else
					(List.map multiplier x :: next, prodNext *. sf);;

let rec generateZerosBelow (mat:matrix) idx row colNum = 
	match mat with
		| [] -> []
		| x :: xs ->	 	
			if(idx >= 0) then 
				x :: generateZerosBelow xs (idx - 1) row colNum
			else if( (List.nth (List.nth mat 0) colNum) = 0.0) then 
				x :: generateZerosBelow xs (idx - 1) row colNum
			else 
				(subtract x row) :: generateZerosBelow xs (idx - 1) row colNum;;

let rec generateGaussForm (mat:matrix) idx n prodAcc= 	
	if(idx < n) then 
		let pvtIndex = findPivotIndex mat 0 idx in
		let pivotRow = List.nth mat pvtIndex in
		let pivot = List.nth (List.nth mat pvtIndex) idx in
		let (swappedMatrix, swapped) = swapRows mat idx pvtIndex in
		let (scaledMatrix, prod) = scaleRowsBelow swappedMatrix 0 idx pivot in
		let zeroedMatrix = generateZerosBelow scaledMatrix idx pivotRow idx in

		generateGaussForm zeroedMatrix (idx + 1) n (prodAcc *. prod *. swapped)
	else
		(mat, prodAcc);;

let rec diagProd (mat:matrix) idx = 
	match mat with
		| [] -> 1.0
		| x :: xs -> (List.nth x idx) *. (diagProd xs (idx + 1));;

let rec detm (m:matrix): float =
	let 
		(gaussForm, prod) = (generateGaussForm m 0 (List.length m) 1.0)
	in
	(diagProd gaussForm 0) /. prod;;

(*
	Inverse calculation also uses Gauss-Jordan elimination. In contrast
	to the previous algorithm, in this, all column are scaled, not just 
	the ones below the diagonal. This leads to the conversion to a
	diagonal matrix. 
	After that, we proceed to scaling the terms of the main diagonal. 
	The same operations are applied on an identity matrix alongside. The
	final state of the identity matrix indicates the inverse.
	
	In this too, complexity is O(n^3) due to the same arguments advanced
	for the determinant.
	Test cases:
		1) [[1.; 3.; 0.]; [2.; 5.; 1.]; [1.; 2.; 1.]]
		2) [[13.; 14.; 15.; 16.]; [1.; 2.; 3.; 4.]; [5.; 6.; 7.; 8.]; [9.; 10.; 11.; 12.]]
		3) scalarmultm 4. (mkunitm 100)
		4) [[0.; 0.; 0.]; [0.; 0.; 0.]; [0.; 0.; 0.]] 
*)

let rec swapRowsIden (mat:matrix) (iden:matrix) r1 r2 = 	
	let row1Mat = List.nth mat r1 in
	let row2Mat = List.nth mat r2 in
	let row1Iden = List.nth iden r1 in
	let row2Iden = List.nth iden r2 in
	if(r1 = r2) then 
		(mat, iden, 1.0)
	else
		(swap mat r1 r2 row1Mat row2Mat, swap iden r1 r2 row1Iden row2Iden, (-1.0));;

let rec multRowOps (mat:matrix) (iden:matrix) idx col pivot = 
	match mat, iden with
		| [], [] -> ([], [], 1.0)
		| [], y :: ys -> raise UnequalMatrixShape
		| x :: xs, [] -> raise UnequalMatrixShape
		| x :: xs, y :: ys -> 	
			if(pivot = 0.0) then 
				(mat, iden, 1.0)
			else
				let denom = List.nth x col in
				let sf = scaleFactor pivot denom in
				let multiplier = mult sf in
				let (next, nextIden, prodNext) = multRowOps xs ys (idx + 1) col pivot in
				if(idx = col) then
					(x :: next, y :: nextIden, prodNext)
				else
					(List.map multiplier x :: next, List.map multiplier y :: nextIden, prodNext *. sf);;


let rec generateZerosInColumn (mat:matrix) (iden:matrix) idx rowMat rowIden colNum = 
	match mat, iden with
		| [], [] -> ([], [])
		| [], y :: ys -> raise UnequalMatrixShape
		| x :: xs, [] -> raise UnequalMatrixShape
		| x :: xs , y :: ys->	
			let (nextMat, nextIden) = generateZerosInColumn xs ys (idx + 1) rowMat rowIden colNum in
			if(idx = colNum) then 
				(x :: nextMat, y :: nextIden)
			else if( (List.nth (List.nth mat 0) colNum) = 0.0) then 	
				(x :: nextMat, y :: nextIden)
			else 
				((subtract x rowMat) :: nextMat, (subtract y rowIden) :: nextIden);;

let rec makeIdentity (matRow:vector) (idenRow:vector) idx = 
	let diagVal = List.nth matRow idx in
	let multiplier = mult (1.0 /. diagVal) in

	if(diagVal = 0.0) then 
		raise SingularMatrix
	else 
		(List.map multiplier matRow, List.map multiplier idenRow);;

let rec makeIdentityMatrix (mat:matrix) (iden:matrix) idx = 
	match mat, iden with
		| [], [] -> ([], [])
		| [], y :: ys -> raise UnequalMatrixShape
		| x :: xs, [] -> raise UnequalMatrixShape
		| x :: xs, y :: ys -> 	
			let (mRow, iRow) = makeIdentity x y idx in
			let (mNext, iNext) = makeIdentityMatrix xs ys (idx + 1)in
			(mRow :: mNext, iRow :: iNext);;

let rec generateInverseIter (mat:matrix) (iden:matrix) idx n prodAcc= 	
		if(idx < n) then 
			let pvtIndex = findPivotIndex mat 0 idx in
			let pivotRowMatrix = List.nth mat pvtIndex in
			let pivotRowIden = List.nth iden pvtIndex in
			let pivot = List.nth (List.nth mat pvtIndex) idx in
			let (swappedMatrix, swappedIden, swapped) = swapRowsIden mat iden idx pvtIndex in
			let (scaledMatrix, scaledIden, prod) = multRowOps swappedMatrix swappedIden 0 idx pivot in
			let (zeroedMatrix, zeroedIden) = generateZerosInColumn scaledMatrix scaledIden 0 pivotRowMatrix pivotRowIden idx in

			generateInverseIter zeroedMatrix zeroedIden (idx + 1) n (prodAcc *. prod *. swapped)
		else
			(mat, iden, prodAcc);;

let rec invm (m:matrix): matrix = 
		let size = List.length m in
		let iden = mkunitm size in
		let (diagForm, diagIden, _) = generateInverseIter m iden 0 size 1.0 in
		let (identity, inverse) = makeIdentityMatrix diagForm diagIden 0 in

		inverse;;





(*
	crossprodv is O(1). This has been implemented considering only three dimensional
	vectors as input.
	Test case:
		1) [] []
		2) [2.; 3.; 4.;] [0.; 1.; -2.]
		3) [2.] [3.; 4.; 4.]
*)
(*let rec crossprodv (v1:vector) (v2:vector): vector = match v1, v2 with
	| [a1; a2; a3], [b1; b2; b3] -> [((a2 *. b3) -. (a3 *. b2)); ((a3 *. b1) -. (a1 *. b3)); ((a1 *. b2) -. (a2 *. b1))] (*standard defn*)
	| x, y -> raise InvalidInput;;*)


let rec subVector i n1 n2 (v:vector) = match v with
		[] 		-> []
	| x :: xs 	->
		if(i >= n1 && i <= n2) then 
			x :: (subVector (i + 1) n1 n2 xs)
		else
			(subVector (i + 1) n1 n2 xs);;

let rec removeColumn (v:vector) i c= 
	(subVector 0 (i + 1) (c - 1) v) @ (subVector 0 0 (i - 1) v);;

let rec cofactor i c (m:matrix) = 
	match m with
			[] 		-> []
		|x :: xs 	-> (removeColumn x i c) :: cofactor i c xs;;

let rec adjoint i c (m:matrix) = 
	if i = c then []
	else 
		if(i mod 2 = 0) then (detm (cofactor i c m)) :: (adjoint (i + 1) c m)
		else -(detm (cofactor i c m)) :: (adjoint (i + 1) c m);;

let rec crossprodv (m:matrix) = 
	let (r, c) = mdim m in 
	adjoint 0 c m;; 


(* 	Code for benchmarking
	Speed test result:
		On an identity matrix of size 500
		Time taken for 
			invm: 8.705s
			detm: 37.605s


	Remarks on the speed test:
		It is not possible to conduct this test on 
		random matrices for 2 reasons:
			1) 	Answers quickly escape to infinity with increasing n
			2) 	Verifying is also difficult due to increasing 
				floating point errors.
		Hence, this has been performed on identity matrices 
		in order to check a minimal condition of efficiency as 
		compared to a regular O(n!) algorithm. 
*) 

let dTime m =
    let t = Sys.time() in
    let fx = detm m in
    Printf.printf "Execution time: %fs\n" (Sys.time() -. t);
    fx


let iTime m =
    let t = Sys.time() in
    let fx = invm m in
    Printf.printf "Execution time: %fs\n" (Sys.time() -. t);
    fx