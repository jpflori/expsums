function FFInv(t)
    if t eq 0 then
        return 0;
    end if;
    return 1/t;
end function;

function KloostermanSum(a)
    K := Parent(a);
    m := Degree(K);
    E := EllipticCurve([1,0,0,0,a]);
    return #E - 2^m;
end function;

function KloostermanSumPartial(a, b)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    pow := ZZ!((2^m-1)/3);
    S := 0;
    t := K!1;
    for i in [0..2^m-2] do
        t *:= z;
        if t^pow eq b then
            S +:= (-1)^(ZZ!(Trace(a*t+FFInv(t))));
        end if;
    end for;
    return S;
end function;

function CubicSum(a)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    S := 0;
    t := K!1;
    for i in [0..2^m-2] do
        t *:= z;
        S +:= (-1)^(ZZ!(Trace(a*t^3)));
    end for;
    return S;
end function;

function CubicSumInv(a, b)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    S := 0;
    t := K!1;
    for i in [0..2^m-2] do
        t *:= z;
        S +:= (-1)^(ZZ!(Trace(a*t^3 + b/t^3)));
    end for;
    return S;
end function;

//this is the mysterious sum, which is up to some squares
// \sum_{u_1 != 1, \psi_3(\tr(u_1)) = b \psi_3(w_2)} \chi(\tr(a w_1 u_1))
//only made for v == 2
function S12U1left(a, b, w1, w2)
    ZZ := Integers();
    K := Parent(w1);
    n := Degree(K);
    m1 := ZZ!(n/2);
    m2 := ZZ!(m1/2);
    z := PrimitiveElement(K);
    //pow := ZZ!((2^n-1)/3);
    pow3 := ZZ!((2^m2+1)/3);
    pow2 := 2^m2;
    pow1 := 2^m1-1;
    t1 := z^pow1;
    t1inv := t1^-1;
    t1pow := t1^-2;
    t1powinv := t1^2;
    u1pow := a*w1^2;
    u1powinv := u1pow;
    u1 := K!w2;
    u1inv := K!w2;
    S := 0;
    //sum for u2 \neq 1
    for i in [0..2^(m1-1)-1] do
        u1pow *:= t1pow;
        u1powinv *:= t1powinv;
        u1 *:= t1;
        u1inv *:= t1inv;
        s2 := u1 + u1inv;
//        s2 := b*(s2^-pow);
        s2 := s2^2;
        s2 := s2^pow3;
        s2pow := s2^pow2;
        if s2pow eq b*s2 then
//        if s2 eq K!1 then
            s1 := (-1)^(ZZ!(Trace(u1pow))) + (-1)^(ZZ!(Trace(u1powinv)));
            S +:= s1;
        end if;
    end for;
    return S;
end function;

// technical stuff to go through cyclotomic classes mod 2^m1
function iatoi(W, n)
    a := 0;
    for i in [1..n] do
        a := ShiftLeft(a,1);
        a +:= W[i];
    end for;
    return a;
end function;

function InitNecklaces(n)
    return [0 : i in [0..n]];
end function;

procedure NextNecklaces(~W, ~a, ~p, n)
    while true do
        p := n;
        while (W[p] eq 1) do
            p := p - 1;
        end while;
        if (p eq 1) then
            W[1] := 1;
            a := iatoi(W, n);
            return;
        end if;
        W[p] := 1;
        for i in [1..(n - p)] do
            W[i + p] := W[i];
        end for;
        if ((n mod p) eq 0) then
            a := iatoi(W, n);
            return;
        end if;
    end while;
end procedure;

// only made for v == 2
// and checks only when K_m1(a) = 1 mod 3
function CheckConjecture(v, mv)
    m2 := 2^(v-2)*mv;
    m1 := 2^(v-1)*mv;
    m0 := 2^v*mv;
    ZZ := Integers();
    K1<z1> := GF(2,m1);
    K0<z0> := GF(2,m0);
    w1 := z0^(2^m1-1);
    w2 := z1^(2^m2-1);
    b := GF(4)!1;
    d := -2^m2;
    pow2m := ShiftLeft(1, m1);
    pow2m1 := pow2m - 1;
    W := InitNecklaces(m1);
    j := 0;
    p := 0;
    //for j in [1..2^m1-1] do
    repeat
        NextNecklaces(~W, ~j, ~p, m1);
        if (j mod 3) ne 0 then
            continue;
        end if;
        a := z1^j;
        KS := KloostermanSum(a);
        KSm := -(2^(m2+1)-KS)/3;
        q := (2^(m2-1) - 1)/3;
        if ((KS - 1) mod 3) eq 0 then
            print j, KS;
            w := Cputime();
            for i in [1..2^(m1-1)] do
                Sab := S12U1left(a, b, w1^i, K0!1);
                D := Sab + KSm + (1 - (-1)^(ZZ!Trace(a*w1^(2*i)))) * q;
                if D ne d and D ne 0 then
                    print mv, j, i, D;
                    return false;
                end if;
            end for;
            print Cputime(w);
        end if;
    //end for;
    until (j eq pow2m1);
    return true;
end function;

v := 2;
mv := 7;
m2 := 2^(v-2)*mv;
m1 := 2^(v-1)*mv;
m0 := 2^v*mv;
ZZ := Integers();
Kvm<zvm> := GF(2,mv*2);
K1<z1> := GF(2,m1);
K0<z0> := GF(2,m0);
w1 := z0^(2^m1-1);
w2 := z1^(2^m2-1);

//b := Random(GF(4));
b := GF(4)!1;

/*
a := Random(K1);
KS := KloostermanSum(a);
while ((KS - 1) mod 3) ne 0 do
    a := Random(K1);
    KS := KloostermanSum(a);
end while;
*/

// Looking for the mysterious function
d := -2^m2;
//for j in [0] do
for j in [1..2^m1-1] do
    a := z1^j;
    if (j mod 3) ne 0 then
        continue;
    end if;
    KS := KloostermanSum(a);
    KS0 := KloostermanSum(K0!a);
    KSz1 := KloostermanSum(z1^(ZZ!(j/3)));
    //if KS eq -8 then
    if ((KS - 1) mod 3) eq 0 then
        Sm := {};
        Sp := {};
        print j, a, KS, KS0, KSz1;
        for i in [1..2^(m1-1)] do
            tst := KloostermanSum(w1^i);
            Sab := S12U1left(a, b, w1^(i), K0!1);
            D := Sab - (2^(m2+1)-KS)/3 + (1 - (-1)^(ZZ!Trace(a*w1^(2*i)))) * (2^(m2-1) - 1)/3;
            print j, a, i, KS, KS0, D, tst;//, Trace(z1^(ZZ!(j/3))*w1^(i)), (w1^i + w1^-i)^(ZZ!((2^m1-1)/3)), KloostermanSum(K1!(1/(w1^i+w1^-i)));
            if D eq d then
                Include(~Sm, tst);
            else
                Include(~Sp, tst);
            end if;
        end for;
        print #Sp, #Sm, #(Sp meet Sm);
        print Sp meet Sm;
    end if;
end for;
