.func TI_S
TI_S:
push {lr}
push {r4-r7}
nop
nop
nop
nop
nop
nop
nop
nop
   movs r3,#0
   mov r4,r0
   movs r6,#8
   rors r4,r6
   eors r3,r4
   mov r4,r0
   movs r6,#16
   rors r4,r6
   eors r3,r4
   mov r4,r0
   movs r6,#24
   rors r4,r6
   eors r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r0
   movs r6,#16
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#24
   rors r4,r6
   mov r5,r1
   and r4,r5
   eor r3,r4
   mov r4,r0
   mov r5,r0
   movs r6,#16
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r0
   movs r6,#24
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#16
   rors r4,r6
   mov r5,r1
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#24
   rors r4,r6
   mov r5,r1
   movs r6,#8
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   mov r5,r0
   movs r6,#24
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#16
   rors r4,r6
   mov r5,r1
   movs r6,#8
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r1
   movs r6,#16
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   mov r5,r1
   movs r6,#16
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r1
   movs r6,#24
   rors r5,r6
   and r4,r5
   eor r3,r4
   mov r4,r0
   mov r5,r1
   movs r6,#24
   rors r5,r6
   and r4,r5
   eor r3,r4
nop
nop
nop
nop
nop
nop
nop
nop
pop {r4-r7}
pop {pc}
.endfunc
