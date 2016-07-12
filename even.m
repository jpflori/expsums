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

function CubicSum(a, b)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    z3 := z^3;
    t := K!1;
    t3 := K!1;
    CS := 1; // zero
    for i in [0..2^m-2] do
        t *:= z;
        t3 *:= z3;
        CS +:= (-1)^(ZZ!(Trace(a*t^3 + b*t)));
    end for;
    return CS;
end function;

function PartialKloostermanSum(a, c)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    zinv := z^(-1);
    t := K!1;
    tinv := K!1;
    powpsi3 := ZZ!((2^m1-1)/3);
    S := 0;
    for i in [0..2^m-2] do
        t *:= z;
        tinv *:= zinv;
        if t^powpsi3 eq c then
            S +:= (-1)^(ZZ!(Trace(a*t + tinv)));
        end if;
    end for;
    return S;
end function;

function PartialSum(a, c)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    t := K!1;
    powpsi3 := ZZ!((2^m1-1)/3);
    S := 0;
    for i in [0..2^m-2] do
        t *:= z;
        if t^powpsi3 eq c then
            S +:= (-1)^(ZZ!(Trace(a*t)));
        end if;
    end for;
    return S;
end function;

function CubicCharacter(a)
    K := Parent(a);
    m := Degree(K);
    alpha := GF(4)!(a^(ZZ!((2^m - 1)/3)));
    return alpha;
end function;


function RandomKS(K)
    a := Random(K);
    while a eq 0 do
        a := Random(K);
    end while;
    KS := KloostermanSum(a);
    while ((KS - 1) mod 3) ne 0 do
        a := Random(K);
        while a eq 0 do
            a := Random(K);
        end while;
        KS := KloostermanSum(a);
    end while;
    return a;
end function;

function RandomB()
    b := Random(GF(4));
    while b eq 0 do
        b := Random(GF(4));
    end while;
    return b;
end function;

//this is the mysterious sum, which is up to some squares
// \sum_{u_1 != 1, \psi_3(\tr(u_1)) = b \psi_3(w_2)} \chi(\tr(a w_1 u_1))
//only made for v == 2
function S2(m2, a, b, w1, w2)
    ZZ := Integers();
    m1 := 2*m2;;
    m0 := 2*m1;
    K1 := GF(2,m1);
    K0 := GF(2,m0);
    z0 := PrimitiveElement(K0);
    pow3 := ZZ!((2^m2+1)/3);
    pow2 := 2^m2;
    pow1 := 2^m1-1;
    t1 := z0^pow1;
    t1inv := t1^-1;
    u1pow := a*w1^2;
    u1powinv := u1pow;
    u1 := K0!1;
    u1inv := K0!1;
    S := 0;
    psi3 := b*w2^(4*pow3);
    // Got through u1 and u1^-1 at once
    for i in [0..2^(m1-1)-1] do
        // Compute a*w1^2*u1^-1 and a*w1^2*u1
        u1pow *:= t1;
        u1powinv *:= t1inv;
        // Compute u1 and u1^-1
        u1 *:= t1;
        u1inv *:= t1inv;
        // Compute trace
        s2 := (u1 + u1inv);
        // Compute trace^(2*(2^m2+1)/3)
        s2 := s2^pow3;
        if s2^pow2 eq psi3*s2 then
//        if s2 eq K0!1 then
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
// and checks only when w1 = w2 = 1
function CheckTheoremB(mv)
    v := 2;
    m2 := 2^(v-2)*mv;
    m1 := 2^(v-1)*mv;
    m0 := 2^v*mv;
    ZZ := Integers();
    K1<z1> := GF(2,m1);
    K0<z0> := GF(2,m0);
    w1 := z0^(2^m1-1);
    w2 := z1^(2^m2-1);
    pow2m := 2^(m2);
    pow2m1 := 2^m1 - 1;
    powpsi3 := ZZ!(pow2m1/3);
    W := InitNecklaces(m1);
    j := 0;
    p := 0;
    //for j in [1..2^m1-1] do
    repeat
        NextNecklaces(~W, ~j, ~p, m1);
        a := z1^j;
        alpha := GF(4)!(a^powpsi3);
        KS := KloostermanSum(a);
        CS := 0;
        if KS mod 3 ne 1 then
            CS := CubicSum(a, a);
        end if;
        print j, pow2m, CS, KS;
        w := Cputime();
        for b in GF(4) do
            if b eq 0 then
                continue;
            end if;
            abne1 := 1;
            if alpha * b eq 1 then
                abne1 := 0;
            end if;
            ainvbne1 := 1;
            if alpha^2 * b eq 1 then
                ainvbne1 := 0;
            end if;
            //print alpha, b, abne1, ainvbne1;
            //print PartialSum(a, b), 1/3*((2-3*abne1)*pow2m-1);
            //print PartialKloostermanSum(a, b), 1/3*((2-3*ainvbne1)*CS+KS-1);
            Tab := ZZ!(((2-3*abne1)*pow2m-(2-3*ainvbne1)*CS-KS)/3);
            Sab := S2(m2, a, b, 1, 1);
            //print Sab, Tab;
            D := Sab - Tab;
            if D ne 0 then
                print mv, j, b, D;
                return false;
            end if;
        end for;
        print Cputime(w);
    //end for;
    until (j eq pow2m1);
    return true;
end function;

// only made for v == 2
// and checks only when K_m1(a) = 1 mod 3 and b = w2 = 1
function CheckConjectureOne(mv)
    v := 2;
    m2 := 2^(v-2)*mv;
    m1 := 2^(v-1)*mv;
    m0 := 2^v*mv;
    ZZ := Integers();
    K1<z1> := GF(2,m1);
    K0<z0> := GF(2,m0);
    w1 := z0^(2^m1-1);
    w2 := z1^(2^m2-1);
    d := -2^m2;
    q := (2^(m2-1) - 1)/3;
    pow2m := 2^m1;
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
        if ((KS - 1) mod 3) eq 0 then
            print j, KS;
            w := Cputime();
            for i in [1..2^(m1-1)] do
                Sab := S2(m2, a, K0!1, w1^i, K0!1);
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

// only made for v == 2
// and checks only when w2 = 1
function StudyConjectureB(mv)
    v := 2;
    m2 := 2^(v-2)*mv;
    m1 := 2^(v-1)*mv;
    m0 := 2^v*mv;
    ZZ := Integers();
    K1<z1> := GF(2,m1);
    K0<z0> := GF(2,m0);
    w1 := z0^(2^m1-1);
    w2 := z1^(2^m2-1);
    pow2m := 2^(m2);
    pow2m1 := 2^m1 - 1;
    powpsi3 := ZZ!(pow2m1/3);
    qchi := (2^(m2-1) - 1)/3;
    qtrace := ZZ!(2*(2^(m2+1) - 1)/3);
    W := InitNecklaces(m1);
    j := 0;
    p := 0;
    //for j in [1..2^m1-1] do
    repeat
        NextNecklaces(~W, ~j, ~p, m1);
        a := z1^j;
        alpha := GF(4)!(a^powpsi3);
        KS := KloostermanSum(a);
        CS := 0;
        if KS mod 3 ne 1 then
            CS := CubicSum(a, a);
        end if;
        print j, pow2m, CS, KS;
        w := Cputime();
        for b in GF(4) do
            //if b eq 0 then
            if b ne 1 then
                continue;
            end if;
            abne1 := 1;
            if alpha * b eq 1 then
                abne1 := 0;
            end if;
            ainvbne1 := 1;
            if alpha^2 * b eq 1 then
                ainvbne1 := 0;
            end if;
            //print alpha, b, abne1, ainvbne1;
            //print PartialSum(a, b), 1/3*((2-3*abne1)*pow2m-1);
            //print PartialKloostermanSum(a, b), 1/3*((2-3*ainvbne1)*CS+KS-1);
            Tab := ZZ!(((2-3*abne1)*pow2m-(2-3*ainvbne1)*CS-KS)/3);
            //Tab := ZZ!(((2)*pow2m-(2)*CS-KS)/3);
            for i in [0..2^(m1-1)] do
                Sab := S2(m2, a, b, w1^i, 1);
                //print Sab, Tab;
                D := Sab - Tab;
                //printf "%3o ", D;
                printf "%3o ", D;
            end for;
            printf "\n";
        end for;
        print Cputime(w);
    //end for;
    until (j eq pow2m1);
    return true;
end function;

// only made for v == 2
// and checks only when K_m1(a) = 1 mod 3
function CountW1(m2, a, b, w2)
    m1 := 2*m2;
    m0 := 2*m1;
    ZZ := Integers();
    K0<z0> := GF(2,m0);
    w1 := z0^(2^m1-1);
    d := -2^m2;
    KS := KloostermanSum(a);
    KSm := -(2^(m2+1)-KS)/3;
    q := (2^(m2-1) - 1)/3;
    if ((KS - 1) mod 3) ne 0 then
        return false;
    end if;
    trz := 0;
    tro := 0;
    for i in [0..2^(m1)] do
        Sab := S2(m2, a, b, w1^i, w2);
        D := Sab + KSm + (1 - (-1)^(ZZ!Trace(a*w1^(2*i)))) * q;
        tr := ZZ ! (D / d);
        if tr eq 0 then
            trz +:= 1;
        else
            if tr ne 1 then
                return false;
            end if;
            tro +:= 1;
        end if;
    end for;
    //c := b*w2^(ZZ!((2^m1-1)/3));
    c := b*w2^(ZZ!(4*(2^m2+1)/3));
    if c eq 1 then
        print c, KS, trz, 2^(m1-1) + 5/6*(KS-4) + 3;
    else
        print c, KS, trz, 2^(m1-1) - 1/6*(KS-4);
    end if;
    return true;
end function;

// only made for v == 2
// and checks only when K_m1(a) = 1 mod 3
function StudyTrace(m2, J0, J1, I0, I1)
    m1 := 2*m2;
    m0 := 2*m1;
    ZZ := Integers();
    K0<z0> := GF(2,m0);
    K1<z1> := GF(2,m1);
    w1 := z0^(2^m1-1);
    d := -2^m2;
    qtrace := ZZ!(2*(2^(m2+1) - 1)/3);
    qchi := 2*(2^(m2-1) - 1)/3;
    psi3 := ZZ!((2^m1-1)/3);
    KK<beta> := GF(4);
    for j in [J0..J1] do
        if j mod 3 ne 0 then
            continue;
        end if;
        a := z1^j;
        KS := KloostermanSum(a);
        if ((KS - 1) mod 3) ne 0 then
            continue;
        end if;
        KSm := -(2^(m2+1)-KS)/3;
        for i in [I0..I1] do
            printf "%2o %5o | %2o |", j, KS, i;
            for k in [0..2] do
                b := beta^k;
                Sab := S2(m2, a, b, w1^i, 1);
                Dchi := Sab + KSm + (ZZ!Trace(a*w1^(2*i))) * qchi;
                Dtrace := Sab + KSm + (ZZ!Trace(a*w1^(2*i))) * qtrace;
                trchi := ZZ ! (Dchi / d);
                //printf " %2o", trchi;
                trtrace := (-1)^(ZZ!Trace(a*w1^(2*i))) *(ZZ ! (Dtrace / d));
                printf " %2o", trtrace;
                //if trtrace eq 0 then
                //    printf "  %1o", k;
                //end if;
            end for;
            printf " #";
            printf " %o", (1/a*K1!(w1^(2*i)+w1^-(2*i)))^psi3;
            //printf " %o", (a*K1!(w1^(2*i)+w1^-(2*i)))^psi3;
            //printf " %o", Trace(a,GF(4));//+GF(4)!((K1!(w1^(2*i)+w1^-(2*i)))^psi3);
            //printf " %o", Trace(z1^(ZZ!(j/3)), GF(4))*(w1^(2*i)+w1^(-2*i))^psi3;
            //printf " %o", Trace(z1^(ZZ!(j/3)), GF(4));//*w1^(2*i)), GF(4));
            //printf " %o", Log(beta, KK!( ((a*(w1^(2*i)+w1^-(2*i))))^psi3));
            //printf " %o", KloostermanSum(a*(w1^(2*i)+w1^-(2*i)));
            printf "\n";
        end for;
    end for;
    return true;
end function;

mv := 3;
v := 2;
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
// a := Random(K1);
a := RandomKS(K1);
